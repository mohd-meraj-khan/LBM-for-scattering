#!/bin/bash

DATA_DIR="data"

#mkdir -p "$DATA_DIR"



a_values=(10 20 30 40 50 60 70 80 90 100)



for ((i = 0; i < ${#a_values[@]}; i++)); do
    a="${a_values[$i]}"
    
	python3 MainRadForce.py $a  # Execute the program with the parameters
done



