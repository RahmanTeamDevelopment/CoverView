#!/bin/bash
unset PYTHONPATH

CFLAGS="-fgnu89-inline"
ABSOLUTE_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PIP_ARGS='--no-cache-dir --ignore-installed --force-reinstall'

virtualenv -p python2.7 --no-site-packages --always-copy env

source ${ABSOLUTE_PATH}/env/bin/activate

pip install --no-cache-dir --ignore-installed --force-reinstall --upgrade pip
pip install ${PIP_ARGS} -r requirements.txt
pip install -U .

check_installation_succeeded.bash
