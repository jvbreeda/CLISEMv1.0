#!/bin/bash

#SCRIPT TO CALCULATE THE JANUARY TEMPERATURE FOR THE MONTE CARLO SIMULATIONS
#THIS SCRIPT MAKES USE OF THE FILE emulice_full_dt500.txt, CONTAINING THE ICE SHEET PARAMATER (ICE VOLUME/ICE AREA) FOR THE ENTIRE LENGTH OF THE EXPERIMENT


index=1
while read line; 
do 

echo $line > EMULICE

echo ${index} x1000 years

Rscript prediction_jan_montecarlo_dt500.R
cp -r tempjan.nc temp_jan_dt500_${index}.nc

let index=index+1

done < emulice_full_dt500.txt

