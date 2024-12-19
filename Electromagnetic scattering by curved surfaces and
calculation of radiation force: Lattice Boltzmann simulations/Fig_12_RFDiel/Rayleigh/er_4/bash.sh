#!/bin/bash

DATA_DIR="data"





ratio_values=(0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)



for ((i = 0; i < ${#ratio_values[@]}; i++)); do
    ratio="${ratio_values[$i]}"
    
	python3 MainRadForce.py $ratio  # Execute the program with the parameters
done


