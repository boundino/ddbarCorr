#!/bin/bash

g++ skim.cc $(root-config --libs --cflags) -I"../includes/" -g -o skim.exe || exit 1

for ii in 18 19 20
do
    inputfile=/eos/user/c/caber/D0D0barAnalysis2018PbPbTrees/d0ana_PbPb2018_HIMinimumBias_${ii}.root
    outputfile=/eos/cms/store/group/phys_heavyions/wangj/DntupleRun2018/skim_d0ana_PbPb2018_HIMinimumBias_${ii}.root
    [[ ${1:-0} -eq 1 ]] && ./skim.exe $inputfile $outputfile
done

rm skim.exe &> /dev/null
