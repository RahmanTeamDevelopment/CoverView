#!/bin/bash
unset PYTHONPATH
CFLAGS="-fgnu89-inline"
ABSOLUTE_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PIP_ARGS='--no-cache-dir --ignore-installed --force-reinstall --no-binary :all: -v -v -v'
virtualenv -p python2.7 --no-site-packages --always-copy env
source ${ABSOLUTE_PATH}/env/bin/activate
pip install --no-cache-dir --ignore-installed --force-reinstall --upgrade pip
source ${ABSOLUTE_PATH}/env/bin/activate

pip install ${PIP_ARGS} Cython==0.25.2
pip install ${PIP_ARGS} six==1.10.0
pip install ${PIP_ARGS} appdirs==1.4.3
pip install ${PIP_ARGS} packaging==16.8
pip install ${PIP_ARGS} py==1.4.32
pip install ${PIP_ARGS} pyparsing==2.2.0
pip install ${PIP_ARGS} pysam==0.10.0
pip install ${PIP_ARGS} pytest==3.0.7
pip install -U .

check_installation_succeeded.bash
