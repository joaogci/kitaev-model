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
module purge
module load gcc/14.2.0
module load gsl
module list

cd src
make 
cd ..

echo "[3] -  Successful compilation"



