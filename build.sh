#!/bin/bash 
PROGNAME=$0

usage() {
    cat << EOF >&2

Usage               : $PROGNAME [option]

option              : nothing for now...

EOF
}

echo "[1] - Exporting directories"
export KITAEV_DIR=$(pwd)

echo "[2] - Compiling SSE code"
cd src
make 
cd ..
echo "[3] -  Successful compilation"



