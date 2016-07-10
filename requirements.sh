#!/bin/bash

set -e

if conda --version >/dev/null 2>&1
then
    conda install astropy mock numpy sphinx pytz lxml html5lib beautifulsoup4
else
    if ! python --version | grep -q 2.6
    then pip install astropy
    fi
    pip install mock numpy sphinx pytz
fi
pip install de405==1997.1 de423==2010.1
pip install https://github.com/brandon-rhodes/assay/archive/master.zip
pip install -e .
