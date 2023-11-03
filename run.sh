#!/bin/sh

date
export LD_LIBRARY_PATH="/home/$USER/Documents/Notes/CellularAutomata_Fast/lib"
make clean
make
bin/out
