#!/bin/bash

g++ skim_mc.cc $(root-config --libs --cflags) -I"../includes/" -g -o skim_mc.exe || exit 1

inputfile=/eos/cms/store/group/phys_heavyions/wangj/DntupleRun2018/mc_analysisTree_promptD0_official.root
outputfile=/eos/cms/store/group/phys_heavyions/wangj/DntupleRun2018/skim_mc_analysisTree_promptD0_official.root
[[ ${1:-0} -eq 1 ]] && ./skim_mc.exe $inputfile $outputfile

rm skim_mc.exe &> /dev/null
