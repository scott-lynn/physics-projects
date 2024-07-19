#!/bin/bash
#Script that allows the user to do a single run with specified parameters
#Asks for user input of params

#Compile program using gfortran with O3 optimisation
gfortran -O3 duffosc.exe -std=f2003 duffosc_v5.f90

#Set A and B constant
A=1.0
B=1.0

#Ask if user wants to run Runge-Kutta testing?
echo 'Run Runge-Kutta testing? y/n'
read switch
if [ $switch == "y" ]
then
    switch=1
    echo 'Testing enabled'
else 
    switch=0
    echo 'Testing disabled'
fi

#Send decision output to file
rm -f switch.dat
echo $switch >> switch.dat

#Ask user for input and read in params
echo Duffing Oscillator
echo Enter inital parameters C, D, Omega as a real number e.g 1.0
echo C:
read C
echo D:
read D
echo omega:
read omega

#Store params in file params.dat
rm -f params.dat 
echo $A > params.dat
echo $B >> params.dat
echo $C >> params.dat
echo $D >> params.dat
echo $omega >> params.dat

./duffosc.exe 