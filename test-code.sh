#!/bin/bash

cd "$(readlink -f $(dirname "${BASH_SOURCE[0]}"))"/ci
echo 'Changing to CI directory: cd' $(pwd)

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
r=$(git rev-parse --show-toplevel)
if grep ' $' \
        $(git ls-files $r/design $r/examples $r/skyfield | grep '\.py$') \
        /dev/null  # prevent hanging on a grep of stdin if ls-files fails
then
    echo
    echo 'Error: trailing whitespace detected on the above-listed lines'
    exit 1
fi
if python --version | grep -q 'Python 3' && command -v pyflakes >/dev/null
then
    d=$(python -c 'import skyfield as s; print(s.__file__.rsplit("/", 1)[0])')
    pyflakes $(find "$d" -name '*.py')
fi
exec assay --batch skyfield.tests
