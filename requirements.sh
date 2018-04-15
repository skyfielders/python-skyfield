#!/bin/bash

set -e

if conda --version >/dev/null 2>&1
then
    conda install \
          astropy=3.0.1 \
          mock=2.0.0 \
          numpy=1.14.2 \
          sphinx=1.7.2 \
          pytz \
          lxml=4.2.1 \
          html5lib=1.0.1 \
          beautifulsoup4=4.6.0 \

else
    if ! python --version 2>&1 | grep -q 2.6
    then pip install astropy
    fi
    pip install mock numpy sphinx pytz
fi
pip install de405==1997.1 de423==2010.1
pip install https://github.com/brandon-rhodes/assay/archive/master.zip
pip install -e .
