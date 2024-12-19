#!/bin/bash

DATA_DIR="data"

#mkdir -p "$DATA_DIR"


# Rayleigh regime
#ratio_values=(0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)

# Mie regime
#ratio_values=(0.5 1.0 1.5 2.0)

# GO regime
ratio_values=(2.5 3.0 3.5 4.0)



for ((i = 0; i < ${#ratio_values[@]}; i++)); do
    ratio="${ratio_values[$i]}"
    
	python3 MainRadForce.py $ratio  # Execute the program with the parameters
done



