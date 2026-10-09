#!/usr/bin/env python3
"""Run post-hoc Experiment 3b against the read-only MRI checkout."""
import sys
from pathlib import Path
from candidate_sensitivity import main

if __name__ == '__main__':
    default_out = Path(__file__).resolve().parents[2] / 'reproduced' / 'reviewer_checks' / 'experiment_3' / 'experiment_3b'
    if '--output-root' not in sys.argv:
        sys.argv += ['--output-root', str(default_out)]
    sys.argv.insert(1, '3b')
    main()
