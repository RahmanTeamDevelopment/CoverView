#!/bin/bash

virtualenv -p python2.7 --no-site-packages --always-copy env
source env/bin/activate
pip install -U pip
pip install -r requirements.txt --no-cache-dir --ignore-installed
pip install . --no-cache-dir --ignore-installed
