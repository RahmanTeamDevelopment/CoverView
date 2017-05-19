#!/bin/bash

virtualenv -p python2.7 --no-site-packages --always-copy env
source env/bin/activate
pip install --no-cache-dir --ignore-installed --force-reinstall --upgrade pip
source env/bin/activate
pip install -r requirements.txt --no-cache-dir --ignore-installed --force-reinstall
pip install . --no-cache-dir --ignore-installed --force-reinstall
