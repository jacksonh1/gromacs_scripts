"""Tests for ExchangeStats — the GROMACS 'Replica exchange statistics' block."""

import pytest

from gromd_analysis.remd_log import ExchangeStats, RemdLogError


def make_log(n_replicas=4, attempts=(999, 500, 499), probs=None, n_exc=None, avg_exc=None):
    """A minimal replica log carrying only what ExchangeStats.parse reads."""
    n_pairs = n_replicas - 1
    probs = probs if probs is not None else [0.30, 0.31, 0.32][:n_pairs]
    n_exc = n_exc if n_exc is not None else [150, 155, 160][:n_pairs]
    avg_exc = avg_exc if avg_exc is not None else [0.30, 0.31, 0.32][:n_pairs]
    total, n_odd, n_even = attempts
    idx = "Repl  " + " ".join(f"{i:4d}" for i in range(n_pairs))

    def section(label, values):
        vals = "Repl   " + " ".join(f"{v}" for v in values)
        return f"Repl  {label}:\n{idx}\n{vals}"

    return "\n".join([
        f"Repl  There are {n_replicas} replicas:",
        "Replica exchange interval: 500",
        "",
        "Replica exchange statistics",
        f"Repl  {total} attempts, {n_odd} odd, {n_even} even",
        section("average probabilities", probs),
        section("number of exchanges", n_exc),
        section("average number of exchanges", avg_exc),
        "",
    ])


def test_parses_header_and_counts():
    s = ExchangeStats.parse(make_log())
    assert s.n_replicas == 4
    assert s.replex_interval == 500
    assert s.total_attempts == 999
    assert s.n_odd == 500
    assert s.n_even == 499


def test_per_pair_lists_are_parallel_and_one_shorter_than_the_replica_count():
    s = ExchangeStats.parse(make_log(n_replicas=4))
    assert s.n_pairs == 3
    assert len(s.avg_prob) == len(s.n_exchanges) == len(s.acceptance_rate) == 3
    assert len(s.attempts_per_pair) == 3


def test_exchange_counts_are_ints_not_floats():
    s = ExchangeStats.parse(make_log())
    assert s.n_exchanges == [150, 155, 160]
    assert all(isinstance(x, int) for x in s.n_exchanges)


def test_attempts_per_pair_follows_step_parity():
    # Even-indexed pairs are attempted on even-type steps, odd on odd-type.
    s = ExchangeStats.parse(make_log(attempts=(999, 500, 499)))
    assert s.attempts_per_pair == [499, 500, 499]


def test_stats_are_read_from_the_block_not_recounted():
    # Acceptance comes straight from 'average number of exchanges'.
    s = ExchangeStats.parse(make_log(avg_exc=[0.11, 0.22, 0.33]))
    assert s.acceptance_rate == [0.11, 0.22, 0.33]


def test_frozen_dataclass_rejects_mutation():
    s = ExchangeStats.parse(make_log())
    with pytest.raises(Exception):
        s.n_replicas = 99


# ── failure modes ────────────────────────────────────────────────────────────

def test_unfinished_run_raises():
    log = "Repl  There are 4 replicas:\nReplica exchange interval: 500\n"
    with pytest.raises(RemdLogError, match="Has the simulation completed"):
        ExchangeStats.parse(log)


def test_missing_replica_count_raises():
    with pytest.raises(RemdLogError, match="replica count"):
        ExchangeStats.parse(make_log().replace("There are 4 replicas", "xx"))


def test_missing_exchange_interval_raises():
    with pytest.raises(RemdLogError, match="exchange interval"):
        ExchangeStats.parse(make_log().replace("Replica exchange interval: 500", "xx"))


def test_missing_section_raises():
    with pytest.raises(RemdLogError, match="number of exchanges"):
        ExchangeStats.parse(make_log().replace("Repl  number of exchanges:", "Repl  gone:"))


def test_pair_count_mismatch_raises():
    # A truncated block would otherwise produce a silently short per-pair table.
    log = make_log(n_replicas=4, probs=[0.3, 0.31], n_exc=[1, 2], avg_exc=[0.3, 0.31])
    with pytest.raises(RemdLogError, match="expected 3"):
        ExchangeStats.parse(log)
