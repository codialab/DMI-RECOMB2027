#!/usr/bin/env python3
"""Run post-hoc Experiment 4 against the read-only MRI checkout."""
import sys
from pathlib import Path
from candidate_sensitivity import main

if __name__ == '__main__':
    default_out = Path(__file__).resolve().parents[2] / 'reproduced' / 'reviewer_checks' / 'experiment_4'
    if '--output-root' not in sys.argv:
        sys.argv += ['--output-root', str(default_out)]
    sys.argv.insert(1, '4')
    main()
