#!/bin/bash
unset PYTHONPATH
ABSOLUTE_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

virtualenv -p python2.7 --no-site-packages --always-copy env
source ${ABSOLUTE_PATH}env/bin/activate
pip install --no-cache-dir --ignore-installed --force-reinstall --upgrade pip
source ${ABSOLUTE_PATH}env/bin/activate
pip install -r requirements.txt --no-cache-dir --ignore-installed --force-reinstall
pip install -U .
check_installation_succeeded.bash
