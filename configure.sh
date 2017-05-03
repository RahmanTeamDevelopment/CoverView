#!/bin/bash

virtualenv -p python2.7 env
source env/bin/activate
pip install -r requirements.txt
pip install -e .
