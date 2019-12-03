#!/bin/bash

# run=(0 1 2 3 4 5)
run=(0 7 8)
##       0   1   2   3   4   5   6   7   8 9 10 11
pt1min=( 5   5   5   5   5   20  3   5   5)
pt1max=( 999 999 999 999 999 999 999 999 999)
pt2min=( 1   1   4   7   15  1   1   1   1)
pt2max=( 999 4   7   15  30  999 999 999 999)
dy=(     1   1   1   1   1   1   1   2   0.5)
cmin=(   0   0   0   0   0   0   0   0   0)
cmax=(   80  80  80  80  80  80  80  80  80)
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

rm ddbar_fithist.exe
rm ddbar_savehist.exe
