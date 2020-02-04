#!/bin/bash

run=(0)
##       0
pt1min=( 4)
pt1max=( 999)
pt2min=( 2)
pt2max=( 999)
dy=(     1)
cmin=(   0)
cmax=(   80)
label=2pi
##
inputdata=/export/d00/scratch/jwang/DntupleRun2018/skim_d0ana_PbPb2018_HIMinimumBias.root
inputmc=rootfiles/masstpl_PbPb.root

RUN_SAVEHIST=${1:-0}
[[ $RUN_SAVEHIST -eq 1 || $# == 0 ]] && { g++ ddbar_savehist.cc $(root-config --libs --cflags) -I"../includes/" -g -o ddbar_savehist.exe || exit 1 ; }
RUN_FITHIST=${2:-0}
[[ $RUN_FITHIST -eq 1 || $# == 0 ]] && { g++ ddbar_fithist.cc $(root-config --libs --cflags) -I"../includes/" -g -o ddbar_fithist.exe || exit 1 ; }

for rr in ${run[@]}
do
    outputdir=dd_cent${cmin[rr]}-${cmax[rr]}_pt${pt1min[rr]}-${pt1max[rr]}-pt${pt2min[rr]}-${pt2max[rr]}_y${dy[rr]/'.'/'p'}_${label}
    [[ $RUN_SAVEHIST -eq 1 ]] && { ./ddbar_savehist.exe $inputdata $inputmc $outputdir ${pt1min[rr]} ${pt1max[rr]} ${pt2min[rr]} ${pt2max[rr]} ${dy[rr]} ${cmin[rr]} ${cmax[rr]} $label ; }
    [[ $RUN_FITHIST -eq 1 ]] && { ./ddbar_fithist.exe "rootfiles/$outputdir/savehist.root" $outputdir ; }
done

rm ddbar_fithist.exe > /dev/null 2>&1
rm ddbar_savehist.exe > /dev/null 2>&1
