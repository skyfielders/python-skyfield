#!/bin/bash
#
# This script should not change directory, because CI strategically runs
# it from somewhere besides the repository root to prevent Skyfield from
# being imported directly from its source tree.

set -e
if ! command -v assay >/dev/null
then
    cat >&2 <<'EOF'
Error: "assay" command not found

Create a virtual environment and run "pip install -r requirements.txt"
to install all of the tools and libraries for Skyfield development.

EOF
    exit 2
fi
if python --version | grep -q 'Python 2.7'
then
    d=$(python -c 'import skyfield as s; print(s.__file__.rsplit("/", 1)[0])')
    pyflakes "$d"/skyfield/*.py
fi
exec assay --batch skyfield/tests
