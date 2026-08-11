# Cell parameter dictionary

`cells.csv` stores stable cell identity and experimental metadata.

`cell_parameters.csv` stores an open-ended, long-form parameter registry. New
parameters may be added without changing its columns. The authoritative model
input remains the versioned TOML configuration and each run's
`resolved_config.json`.

One cell is one scientific parameter combination. One `(cell_id, replicate_id)`
pair is one stochastic run. Reusing replicate IDs across cells creates matched
seed blocks. Never reuse a cell ID after changing a scientific parameter.
