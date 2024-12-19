#!/bin/bash

DATA_DIR="data"

#mkdir -p "$DATA_DIR"


# Rayleigh regime
#ratio_values=(0.04 0.1 0.16)


# Mie regime
ratio_values=(0.22 0.28 0.34 0.4 0.46 0.52 0.58 0.64 0.7 0.76)



for ((i = 0; i < ${#ratio_values[@]}; i++)); do
    ratio="${ratio_values[$i]}"
    
	python3 MainRadForce.py $ratio  # Execute the program with the parameters
done



