#!/bin/bash

DATA_DIR="data"

#mkdir -p "$DATA_DIR"



ratio_values=(0.5 0.95 1.0)



for ((i = 0; i < ${#ratio_values[@]}; i++)); do
    ratio="${ratio_values[$i]}"
    
	python3 Main_Tot_2D.py $ratio  # Execute the program with the parameters
done



