#!/bin/bash

cd /localdisk/home/azvoleff/teamlucc

git pull
sed 's/dplyr/plyr/' <DESCRIPTION >DESCRIPTION
sed '/dplyr/d' <NAMESPACE >NAMESPACE
