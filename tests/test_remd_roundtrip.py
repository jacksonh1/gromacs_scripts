"""Tests for remd_roundtrip — reconstruction of the replica permutation walk.

The critical property is the self-check: per-pair swap counts recovered from the
`Repl ex` lines must equal the log's `number of exchanges` block. These tests
synthesize a full log (exchange records + a matching statistics block) so the
reconstruction, round-trip counting, dwell, and both self-check outcomes are
exercised without reading any real trajectory.
"""

import pytest

from gromd_analysis.remd_roundtrip import compute, _parse_replex_line
from gromd_analysis.remd_log import RemdLogError


def make_log(n_replicas, attempts, n_exc=None, tamper_exc=None):
    """Full replica log for `attempts` — a list where each element is the list of
    lower-slot indices that swapped on that exchange attempt.

    The `number of exchanges` block is computed from `attempts` so the log is
    self-consistent, unless `tamper_exc` overrides it (to drive the self-check
    failure path).
    """
    n_pairs = n_replicas - 1
    true_exc = [0] * n_pairs
    for swaps in attempts:
        for lo in swaps:
            true_exc[lo] += 1

    lines = [
        f"Repl  There are {n_replicas} replicas:",
        "Replica exchange interval: 500",
        "",
    ]
    for step, swaps in enumerate(attempts, start=1):
        lines.append(f"Replica exchange at step {step*500} time {step*1.0:.5f}")
        cols = []
        for slot in range(n_replicas):
            cols.append(f"{slot:4d}")
            if slot in swaps:            # x sits between slot and slot+1
                cols.append("x")
        lines.append("Repl ex " + " ".join(cols))

    idx = "Repl  " + " ".join(f"{i:4d}" for i in range(n_pairs))
    block_exc = tamper_exc if tamper_exc is not None else (n_exc or true_exc)
    probs = [0.3] * n_pairs
    avg_exc = [0.3] * n_pairs

    def section(label, values):
        return f"Repl  {label}:\n{idx}\nRepl   " + " ".join(str(v) for v in values)

    total = len(attempts)
    lines += [
        "Replica exchange statistics",
        f"Repl  {total} attempts, {total} odd, {total} even",
        section("average probabilities", probs),
        section("number of exchanges", block_exc),
        section("average number of exchanges", avg_exc),
        "",
    ]
    return "\n".join(lines)


def write_log(tmp_path, text):
    log = tmp_path / "rest2.log"
    log.write_text(text)
    return log


def test_reconstructs_one_round_trip(tmp_path):
    # n=3; drive walker 0 through slot 0->1->2->1->0 = exactly one round trip.
    attempts = [[0], [1], [1], [0]]
    log = write_log(tmp_path, make_log(3, attempts))
    mix = compute(log)
    assert mix.n_replicas == 3
    assert mix.n_attempts == 4
    assert mix.round_trips[0] == 1


def test_partial_excursion_is_not_a_round_trip(tmp_path):
    # Walker 0 goes 0->1->0 (never reaches the top slot 2): no round trip.
    attempts = [[0], [0]]
    log = write_log(tmp_path, make_log(3, attempts))
    mix = compute(log)
    assert mix.round_trips[0] == 0


def test_dwell_fractions_sum_to_one(tmp_path):
    attempts = [[0], [1], [1], [0]]
    log = write_log(tmp_path, make_log(3, attempts))
    mix = compute(log)
    for row in mix.dwell_frac:
        assert abs(sum(row) - 1.0) < 1e-9


def test_run_ns_from_last_exchange_time(tmp_path):
    attempts = [[0], [1], [1], [0]]  # last time = 4.0 ps
    log = write_log(tmp_path, make_log(3, attempts))
    mix = compute(log)
    assert mix.run_ns == pytest.approx(0.004)


def test_self_check_passes_on_consistent_log(tmp_path):
    # Larger walk; block counts derived from the same attempts must agree.
    attempts = [[0, 2], [1], [0, 2], [1], [2]]
    log = write_log(tmp_path, make_log(4, attempts))
    mix = compute(log)  # would raise if the self-check failed
    assert mix.n_replicas == 4


def test_self_check_fails_when_block_disagrees(tmp_path):
    # Repl ex lines give pair counts [2,2]; tamper the block to [2,3].
    attempts = [[0], [1], [1], [0]]
    log = write_log(tmp_path, make_log(3, attempts, tamper_exc=[2, 3]))
    with pytest.raises(RemdLogError, match="self-check FAILED"):
        compute(log)


def test_no_replex_lines_raises(tmp_path):
    log = tmp_path / "rest2.log"
    log.write_text("Repl  There are 3 replicas:\nReplica exchange interval: 500\n")
    with pytest.raises(RemdLogError, match="no 'Repl ex' lines"):
        compute(log)


def test_parse_replex_line_wrong_column_count_raises():
    # Line lists 2 slot columns but 3 replicas expected.
    with pytest.raises(RemdLogError, match="expected 3"):
        _parse_replex_line("Repl ex    0    1", 3)


def test_parse_replex_line_recovers_swap_indices():
    assert _parse_replex_line("Repl ex    0 x    1    2 x    3", 4) == [0, 2]
    assert _parse_replex_line("Repl ex    0    1    2    3", 4) == []
