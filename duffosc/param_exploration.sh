#!/bin/bash
#Script that automates parameter exploration of D, C and Omega
version=1

#Disable testing phase
switch=0
rm -f switch.dat
echo $switch >> switch.dat

#Compile program using gfortran with O3 optimisation
gfortran -O3 duffosc.exe -std=f2003 duffosc_v5.f90

#Driving force (D) exploration
version=1

A=1.0
B=1.0
C=0.1
omega=2.0

#Iterate through D parameter
for D in $(seq 1.0 0.01 2.0)
do
            echo "Cycle num:"${version}

            rm -f params.dat 
            echo $A > params.dat
            echo $B >> params.dat
            echo $C >> params.dat
            echo $D >> params.dat
            echo $omega >> params.dat

            ./duffosc.exe 
            cat poincare.dat > data/poincare_D${version}.dat
            cat params.dat > data/params_D${version}.dat

            #Output as graph to xmgrace
            #xmgrace data/poincare_D${version}.dat -param settings2.par &

            version=$(($version + 1 ))      
done

#Damping force (C) exploration
version=1

A=1.0
B=1.0
D=1.0
omega=2.0

#Iterate through C parameter
for C in $(seq 0.1 0.01 0.3)
do
            echo "Cycle num:"${version}

            rm -f params.dat 
            echo $A > params.dat
            echo $B >> params.dat
            echo $C >> params.dat
            echo $D >> params.dat
            echo $omega >> params.dat

            ./duffosc.exe 
            cat poincare.dat > data/poincare_C${version}.dat
            cat params.dat > data/params_C${version}.dat

            #Output as graph to xmgrace
            #xmgrace data/poincare_C${version}.dat -param settings2.par &

            version=$(($version + 1 ))      
done

#Damping force (C) exploration
version=1

A=1.0
B=1.0
C=0.1
D=1.0

#Iterate through C parameter
for omega in $(seq 2.0 0.01 3.0)
do
            echo "Cycle num:"${version}

            rm -f params.dat 
            echo $A > params.dat
            echo $B >> params.dat
            echo $C >> params.dat
            echo $D >> params.dat
            echo $omega >> params.dat

            ./duffosc.exe 
            cat poincare.dat > data/poincare_w${version}.dat
            cat params.dat > data/params_w${version}.dat

            #Output as graph to xmgrace
            #xmgrace data/poincare_w${version}.dat -param settings2.par &

            version=$(($version + 1 ))      
done

#Notify user that process is complete
echo "SCRIPT FINISH"