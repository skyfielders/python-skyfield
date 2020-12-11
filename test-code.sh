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

echo "assay command found"

current_dir=`pwd`
git_root=`git rev-parse --show-toplevel`

cd $git_root
if grep ' $' $(git ls-files design examples skyfield | grep '\.py$')
then
    echo
    echo 'Error: trailing whitespace detected on the above-listed lines'
    exit 1
fi
cd $current_dir

echo "searching for files with a whitespace done."

if python --version | grep -q 'Python 3.6' && command -v pyflakes >/dev/null
then
    d=$(python -c 'import skyfield as s; print(s.__file__.rsplit("/", 1)[0])')
    pyflakes $(find "$d" -name '*.py')
fi

echo "pyflakes executed"

exec assay --batch skyfield.tests
