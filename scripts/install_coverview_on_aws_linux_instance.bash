#!/bin/bash
sudo yum -y install git
sudo yum -y install gcc
sudo yum -y install zlib-devel
mkdir test
cd test
git clone  https://github.com/RahmanTeamDevelopment/CoverView
cd CoverView
./install.bash