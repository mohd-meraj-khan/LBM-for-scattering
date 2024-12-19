#!/bin/bash

DATA_DIR="data"

#mkdir -p "$DATA_DIR"



#a_values=(25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175)



#for ((i = 0; i < ${#a_values[@]}; i++)); do
#    ratio="${a_values[$i]}"
    
#	python3 MainRadForce.py $ratio  # Execute the program with the parameters
#done



start=25
end=175
increment=5

for ((ratio = start; ratio <= end; ratio += increment)); do
    python3 MainRadForce.py "$ratio" 
done

