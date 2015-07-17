#!/bin/bash

for (( i=0; i<2; i++ ))
do
for (( j=0; j<5; j++ ))
do
for (( k=0; k<6; k++ ))
do
#  echo "NReco_Arm${i}Plane${j}CSC${k}->Fill(NHITS[$i][$j][$k]);"

#echo "RecoX_Arm${i}Plane${j}CSC${k} = std::auto_ptr<TH1D>(new TH1D("RecoXArm${i}Plane${j}CSC${k}","RecoXArm${i}Plane${j}CSC${k}",300,-1500,1500));"
#echo "RecoX_Arm${i}Plane${j}CSC${k}->SetDirectory(0);"
#echo "RecoY_Arm${i}Plane${j}CSC${k} = std::auto_ptr<TH1D>(new TH1D("RecoYArm${i}Plane${j}CSC${k}","RecoYArm${i}Plane${j}CSC${k}",300,-1500,1500));"
#echo "RecoY_Arm${i}Plane${j}CSC${k}->SetDirectory(0);"
#echo "NReco_Arm${i}Plane${j}CSC${k} = std::auto_ptr<TH1D>(new TH1D("NRecoXArm${i}Plane${j}CSC${k}","NRecoXArm${i}Plane${j}CSC${k}",100,-0.5,99.5));"
#echo "NReco_Arm${i}Plane${j}CSC${k}->SetDirectory(0);"

#echo "RecoX_Arm${i}Plane${j}CSC${k}->Write();"
#echo "RecoY_Arm${i}Plane${j}CSC${k}->Write();"
#echo "NReco_Arm${i}Plane${j}CSC${k}->Write();"

echo "std::auto_ptr<TH1D> RecoX_Arm${i}Plane${j}CSC${k};"
echo "std::auto_ptr<TH1D> RecoY_Arm${i}Plane${j}CSC${k};"
echo "std::auto_ptr<TH1D> NReco_Arm${i}Plane${j}CSC${k};"
done
done
done
