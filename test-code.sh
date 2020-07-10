#!/bin/bash
#
# This script should not change directory, because CI strategically runs
# it from somewhere besides the repository root to prevent Skyfield from
# being imported directly from its source tree.

set -e
if python --version | grep -q 'Python 3.6' && command -v pyflakes >/dev/null
then
    d=$(python -c 'import skyfield as s; print(s.__file__.rsplit("/", 1)[0])')
    pyflakes $(find "$d" -name '*.py')
fi
pytest -v --pyargs skyfield
