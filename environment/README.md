# Software environment

The selected plotting cells need NumPy, pandas, SciPy, PyArrow and Matplotlib.
They read existing compressed tables and do not require COBRApy or Gurobi.

From the repository root, create an environment **inside this repository**:

```bash
python -m venv .venv
.venv/bin/python -m pip install -r environment/requirements-figures.txt
.venv/bin/python reproduce.py status
.venv/bin/python reproduce.py checksums
.venv/bin/python reproduce.py validate
.venv/bin/python -m pip install -r environment/requirements-validation.txt
.venv/bin/python -m pytest -p no:cacheprovider -q
```

Python 3.12.3 and the pinned plotting versions were observed in the validation
environment and are recorded in `manifests/provenance/VALIDATION.json`. They do
not establish the versions that originally generated scientific outputs.
Rendering can differ with fonts or Matplotlib versions; numerical validation
uses frozen population and endpoint checks, and input identity uses SHA-256.
Figure 1 also requires Inkscape. Optional upstream reconstruction dependencies
are separate and do not make the absent source data available.

`requirements-upstream.txt` records optional modeling dependencies. The source
Stage-11 file pins Troppo **0.0.7**, despite an audit's abbreviated “0.7” wording;
the actual requirement file is preserved as provenance. Original solver version
claims are limited to that file and source manifests. No solver was imported or
run during migration; a successful installation is not scientific qualification.

The reviewer workflows do not run Gurobi. Repeating the upstream solver analyses
is unsupported because required biological inputs, models, and run artifacts are
not included; do not treat the optional dependency list as a turnkey workflow.

`requirements-validation.txt` adds dependencies for the retained synthetic-core
checks. No package installation is performed as part of migration.
