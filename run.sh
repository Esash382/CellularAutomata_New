#!/bin/sh

date
export LD_LIBRARY_PATH="/home/ashraya/Documents/Notes/CellularAutomata_Fast/lib"
make clean
make
bin/out
