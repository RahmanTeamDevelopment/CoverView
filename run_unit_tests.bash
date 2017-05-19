#!/bin/bash
source env/bin/activate
unset PYTHONPATH
pytest test
