#! /bin/bash

run=(0)

pt1min=( 2 4 6 8 4 6 8 6 8 8)
pt1max=( 999 999 999 999 999 999 999 999 999 999)
pt2min=( 2 2 2 2 4 4 4 6 6 8)
pt2max=( 999 999 999 999 999 999 999 999 999 999)
dy=(     1 1 1 1 1 1 1 1 1 1)
cmin=(   0 0 0 0 0 0 0 0 0 0)
cmax=(   80 80 80 80 80 80 80 80 80 80)
label=4pi

inputdata=DntupleRun2018/skim_d0ana_PbPb2018_HIMinimumBias.root
inputmc=rootfiles/masstpl_PbPb.root
genmc=DntupleRun2018/skim_mc_analysisTree_promptD0_official.root
swapmc=DntupleRun2018/mc_analysisTree_promptD0_official.root
eff=efficiency.root

dpairtree=dptree.root
swaptree=swaptree.root

mkdir -p build
cd build
cmake ..
make || exit
cd ..

RUN_SAVETREE=${1:-0}
RUN_2DFIT=${2:-0}

for rr in ${run[@]}
do
    outputdir=dd_cent${cmin[rr]}-${cmax[rr]}_pt${pt1min[rr]}-${pt1max[rr]}_y${dy[rr]/'.'/'p'}_${label}
    [[ $RUN_SAVETREE -eq 1 ]] && { build/ddbar_savetree $inputdata $inputmc $outputdir ${pt1min[rr]} ${pt1max[rr]} ${pt2min[rr]} ${pt2max[rr]} ${dy[rr]} ${cmin[rr]} ${cmax[rr]} $label $eff $swapmc $mctree; }
    [[ $RUN_2DFIT -eq 1 ]] && { build/ddbar_2dfit $"rootfiles/$outputdir/$dpairtree" "rootfiles/$outputdir/${swaptree}" $outputdir ; }
done
