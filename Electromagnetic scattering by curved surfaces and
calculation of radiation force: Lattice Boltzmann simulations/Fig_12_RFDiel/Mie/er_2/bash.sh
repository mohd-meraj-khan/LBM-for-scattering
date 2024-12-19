#!/bin/bash

DATA_DIR="data"





ratio_values=(0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.00 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.10)



for ((i = 0; i < ${#ratio_values[@]}; i++)); do
    ratio="${ratio_values[$i]}"
    
	python3 MainRadForce.py $ratio  # Execute the program with the parameters
done


