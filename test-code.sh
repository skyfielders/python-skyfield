#!/bin/bash

cd "$(dirname "${BASH_SOURCE[0]}")"

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
then pyflakes skyfield/*.py
fi
exec assay --batch skyfield/tests
