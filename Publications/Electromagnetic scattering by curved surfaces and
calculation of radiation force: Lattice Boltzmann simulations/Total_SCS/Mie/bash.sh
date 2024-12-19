#!/bin/bash

DATA_DIR="data"





ratio_values=(0.90 0.925 0.95 0.975 1.0 1.025 1.05 1.075 1.10)



for ((i = 0; i < ${#ratio_values[@]}; i++)); do
    ratio="${ratio_values[$i]}"
    
	python3 MainScatWidth.py $ratio  # Execute the program with the parameters
done


