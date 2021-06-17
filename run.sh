#!/bin/bash

rm -f energy.log traj.pdb
make clean
make -j
./main
