#!/bin/bash

DATA_DIR="data"



start=25
end=175
increment=5

for ((ratio = start; ratio <= end; ratio += increment)); do
    python3 MainRadForce.py "$ratio" 
done

