#!/bin/bash

#SCRIPT TO CALCULATE THE MONTHLY TEMPERATURE AND PRECIPITATION DURING THE ENTIRE LENGTH OF THE SIMULATIONS
#THIS SCRIPT MAKES USE OF THE FILE emulice_full.txt, CONTAINING THE ICE SHEET PARAMATER (ICE VOLUME/ICE AREA) FOR THE ENTIRE LENGTH OF THE EXPERIMENT

#choose timestep equal to 500, 1000 or 2000 (years)
timestep=1000


index=1

while read line; 
do 

echo $line > EMULICE

echo ${index} x $timestep years

Rscript prediction.R
cp -r MONTHLY_TEMP.nc monthly_temp_dt${timestep}_${index}.nc
cp -r MONTHLY_PRECIP.nc monthly_precip_dt${timestep}_${index}.nc

let index=index+1

done < emulice_full_dt$timestep.txt

