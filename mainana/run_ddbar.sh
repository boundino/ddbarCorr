#!/bin/bash

run=(0)
##       0
pt1min=( 2 4 6 8 4 6 8 6 8 8)
pt1max=( 999 999 999 999 999 999 999 999 999 999)
pt2min=( 2 2 2 2 4 4 4 6 6 8)
pt2max=( 999 999 999 999 999 999 999 999 999 999)
dy=(     1 1 1 1 1 1 1 1 1 1)
cmin=(   0 0 0 0 0 0 0 0 0 0)
cmax=(   80 80 80 80 80 80 80 80 80 80)
label=2pi
##
inputdata=DntupleRun2018/skim_d0ana_PbPb2018_HIMinimumBias.root
inputmc=rootfiles/masstpl_PbPb.root
genmc=DntupleRun2018/skim_mc_analysisTree_promptD0_official.root
eff=efficiency.root

[[ -d DntupleRun2018 ]] || ln -s /data/wangj/DntupleRun2018

run_efficiency=1

while [ "$1" != "" ]; do
    case $1 in
        -e | --efficiency )           shift
                                      run_efficiency=1
                                      ;;
        *)                            break
                                      ;;
    esac
done

if [ "$run_efficiency" = "1" ]; then
    echo "Getting efficiency from MC samples ..."
    g++ efficiency.cc -g -Wall $(root-config --libs --cflags) -I"../includes/" -g -o eff
    ./eff $genmc $eff
    rm eff > /dev/null 2>&1
fi

RUN_SAVEHIST=${1:-0}
[[ $RUN_SAVEHIST -eq 1 || $# == 0 ]] && { g++ ddbar_savehist.cc $(root-config --libs --cflags) -I"../includes/" -g -o ddbar_savehist.exe || exit 1 ; }
RUN_FITHIST=${2:-0}
[[ $RUN_FITHIST -eq 1 || $# == 0 ]] && { g++ ddbar_fithist.cc $(root-config --libs --cflags) -I"../includes/" -g -o ddbar_fithist.exe || exit 1 ; }

for rr in ${run[@]}
do
    outputdir=dd_cent${cmin[rr]}-${cmax[rr]}_pt${pt1min[rr]}-${pt1max[rr]}-pt${pt2min[rr]}-${pt2max[rr]}_y${dy[rr]/'.'/'p'}_${label}
    [[ $RUN_SAVEHIST -eq 1 ]] && { ./ddbar_savehist.exe $inputdata $inputmc $outputdir ${pt1min[rr]} ${pt1max[rr]} ${pt2min[rr]} ${pt2max[rr]} ${dy[rr]} ${cmin[rr]} ${cmax[rr]} $label $eff; }
    [[ $RUN_FITHIST -eq 1 ]] && { ./ddbar_fithist.exe "rootfiles/$outputdir/savehist.root" $outputdir ; }
done

rm ddbar_fithist.exe > /dev/null 2>&1
rm ddbar_savehist.exe > /dev/null 2>&1
