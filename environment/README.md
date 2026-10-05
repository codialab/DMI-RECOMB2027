# Software environment

The selected plotting cells need NumPy, pandas, SciPy, PyArrow and Matplotlib.
They read existing compressed tables and do not require COBRApy or Gurobi.

From the repository root, create an environment **inside this repository**:

```bash
python -m venv .venv
.venv/bin/python -m pip install -r environment/requirements-figures.txt
.venv/bin/python reproduce.py status
```

Pinned plotting versions were observed in the migration environment and are
recorded in `manifests/provenance/VALIDATION.json`. They do not establish the
versions that originally generated scientific outputs. Rendering may differ
with fonts or Matplotlib versions; scientific-table identities use SHA-256.

`requirements-upstream.txt` records optional modeling dependencies. The source
Stage-11 file pins Troppo **0.0.7**, despite an audit's abbreviated “0.7” wording;
the actual requirement file is preserved as provenance. Original solver version
claims are limited to that file and source manifests. No solver was imported or
run during migration; a successful installation is not scientific qualification.

For an explicitly authorized future Gurobi solve on this host, use the existing
`/home/pty/my_env/bin/python`, first check all data/output paths, and request only
the necessary native-host command if the known sandbox HostID failure occurs.
Do not change solver, environment, tolerances or model assumptions to work around it.

`requirements-validation.txt` adds dependencies for the retained synthetic-core
checks. No package installation is performed as part of migration.
