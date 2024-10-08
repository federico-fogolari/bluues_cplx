#!/bin/bash

#Path of the directory in which run_nanoshaper is located
CMD=$(echo "readlink ${BASH_SOURCE[0]}")
if [ ! -z "$CMD" ];  then
SCRIPT_DIR=$(dirname $(${CMD}))
echo "Set LD_LIBRARY_PATH:$SCRIPT_DIR"
export LD_LIBRARY_PATH=$SCRIPT_DIR:$LD_LIBRARY_PATH
fi

CONF_PATH=${1?Error: Surface Configuration not given}

echo "Surface Configuration: $1"
$SCRIPT_DIR/NanoShaper.bin $1

