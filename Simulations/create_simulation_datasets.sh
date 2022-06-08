#!/bin/bash

echo "Choose a condition: "
echo "1. A dynamic group of cells following a bifurcation trajectory"
echo "2. A non-dynamic group of cells"
echo "3. A dynamic group of cells following a bifurcation trajectory and a non-dynamic group of cells"
echo "4. A dynamic group of cells following a bifurcation trajectory and 3 non-dynamic groups of cells"

read option

if [[ $option -gt 0 ]] && [[ $option -lt 5 ]] ; then
    echo "You chose $option"
    echo "Choose number of simulations"
    read n
    echo "Choose capture rate"
    read capture_rate
    Rscript ./code/create_dataset.r $option $n $capture_rate
else
    echo "You chose $option. Please choose a number between 1 and 4"
    exit
fi