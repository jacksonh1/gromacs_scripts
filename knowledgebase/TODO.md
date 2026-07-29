# TODO list

- [pending] Restructure output-file behavior so logs/metadata aren't orphaned from the trajectory data. Right now trajectories live on scratch while the logs and other outputs sit in the declared output folder, so the two are disconnected (and scratch purge leaves dangling symlinks). Make scratch the primary location for the complete, self-contained run — logs included, with the actual trajectories in a `trajectories/` subfolder — then in the declared output folder split copy-vs-symlink by size/durability-value, not by a fixed folder list: **copy** everything small enough to survive a scratch purge cheaply (`analysis/`, final PDB, `parameters.txt`, **the logs**, `.mdp`s, final `.tpr`/`.gro`) and **symlink** only the genuinely large files (trajectories, large intermediates). Rationale: logs are the only record of what a run did once trajectories are purged and cost only a few MB, so they must be on the durable (copied) side — symlinking them makes the cheapest, most-worth-keeping artifact the least durable. Dangling symlinks for the large files on purge are the accepted tradeoff of keeping bulk data on ephemeral storage.

- [pending] Package the Python analysis scripts into a simple installable package (e.g. `pyproject.toml` with a console entry point or importable module), so they can be `pip install`-ed rather than run as loose scripts.

- [pending] (optional) Switch the analysis code from MDAnalysis to mdtraj.









## for claude
- The user will ask you to add an item to the TODO list or ask you what is on the TODO list.
- You will update the TODO list accordingly and provide the current list when requested.
- You can track the status of each item as "pending", "in progress", or "completed".
- update the conventions for this TODO list here as we decide on them.
- format - todo list items should be 1-2 sentences long, and should be clear and concise.
