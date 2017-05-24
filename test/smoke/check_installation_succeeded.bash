#!/bin/bash

echo ''
echo 'Checking CoverView installation'

# Test running from the directory in which the script was invoked
coverview --help &> /dev/null

if [ $? != 0 ]
then
    echo 'Error: CoverView installation failed. Cannot run coverview from directory ' $(pwd)
    exit 1
fi

# Test running from directory that contains the script
ABSOLUTE_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd $ABSOLUTE_PATH

coverview --help &> /dev/null

if [ $? != 0 ]
then
    echo 'Error: CoverView installation failed. Cannot run coverview from directory ' $ABSOLUTE_PATH
    exit 1
fi


echo 'CoverView installation succeeded'
echo ''
