#!/bin/bash

# just add the path to the source code into MATLABPATH
REPOPATH=`pwd`
MCMCPATH=${REPOPATH}/src
EXAMPLEPATH=${REPOPATH}/Examples

export MATLABPATH=${MATLABPATH}:${MCMCPATH}:${EXAMPLEPATH}

