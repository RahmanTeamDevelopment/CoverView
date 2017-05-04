#!/bin/bash

virtualenv -p python2.7 env
source env/bin/activate
pip install -r requirements.txt --no-cache-dir --ignore-installed
pip install -e . --no-cache-dir
