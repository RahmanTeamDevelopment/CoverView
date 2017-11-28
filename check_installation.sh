#!/bin/bash

echo ''
echo 'Checking CoverView installation'

./coverview --help &> /dev/null
./gui --help &> /dev/null
./ensembl_db --help &> /dev/null

if [ $? != 0 ]
then
    echo 'Error: CoverView installation failed.'
    exit 1
fi

echo 'CoverView installation succeeded'
echo ''
