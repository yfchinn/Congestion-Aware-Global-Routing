#!/bin/bash

echo "===== Compiling... ====="
make || { echo "Make failed!"; exit 1; }

CASE_DIR="cases"
EVAL_SCRIPT="$CASE_DIR/eval2008.pl"

inputs=("3d.txt" "adaptec1.gr" "adaptec2.gr" "adaptec3.gr" "adaptec4.gr" "adaptec5.gr" "bigblue1.gr" "bigblue2.gr" "bigblue3.gr" "newblue1.gr" "newblue2.gr" "newblue5.gr" "newblue6.gr")
outputs=("3ds1.txt" "adaptec1_sol.gr" "adaptec2_sol.gr" "adaptec3_sol.gr" "adaptec4_sol.gr" "adaptec5_sol.gr" "bigblue1_sol.gr" "bigblue2_sol.gr" "bigblue3_sol.gr" "newblue1_sol.gr" "newblue2_sol.gr" "newblue5_sol.gr" "newblue6_sol.gr")

for i in "${!inputs[@]}"; do
    IN="$CASE_DIR/${inputs[$i]}"
    OUT="$CASE_DIR/${outputs[$i]}"

    echo -e "\n------------------------------------------"
    echo "Processing Case: ${inputs[$i]}"
    
    ./router "$IN" "$OUT"

    if [ -f "$OUT" ]; then
        echo "Running Evaluation..."
        perl "$EVAL_SCRIPT" "$IN" "$OUT"
    else
        echo "Error: Output file $OUT was not generated."
    fi
done

echo -e "\nAll test cases completed."