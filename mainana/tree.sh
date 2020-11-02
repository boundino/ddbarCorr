#! /bin/bash

run=(0)

pt1min=( 2 4 6 8 4 6 8 6 8 8)
pt1max=( 999 999 999 999 999 999 999 999 999 999)
pt2min=( 2 2 2 2 4 4 4 6 6 8)
pt2max=( 999 999 999 999 999 999 999 999 999 999)
dy=(     1 1 1 1 1 1 1 1 1 1)
cmin=(   0 0 0 0 0 0 0 0 0 0)
cmax=(   80 80 80 80 80 80 80 80 80 80)
label=2pi

RUN_SAVETREE=${1:-0}
RUN_2DFIT=${2:-0}

system=pp
inputdata=/data2/mjpeters/HIZeroBiasAll_skim_newNtuple/skim_merged.root
inputmc=rootfiles/masstpl_PbPb.root
inputmixdata=/data2/mjpeters/HIZeroBiasAll_skim_newNtuple/skim_merged.root
genmc=DntupleRun2018/DDbar/skim_mc_analysisTree_promptD0_official.root
swapmc=DntupleRun2018/DDbar/mc_analysisTree_promptD0_official.root
eff=rootfiles/efficiency.root

dpairtree=dptree.root
swaptree=swaptree.root

[[ -d DntupleRun2018 ]] || ln -s /data/wangj/DntupleRun2018

if [[ "$system" == "PbPb" ]]; then
  sed -i -e "s/SYSTEM_PP/SYSTEM_PBPB/g" CMakeLists.txt
elif [[ "$system" == "pp" ]]; then
  sed -i -e "s/SYSTEM_PBPB/SYSTEM_PP/g" CMakeLists.txt
fi
mkdir -p build
cd build
cmake ..
make || exit
cd ..

[[ -f $eff ]] || build/ddbar_efficiency $genmc $eff

for rr in ${run[@]}
do
    outputdir=$system/dd_cent${cmin[rr]}-${cmax[rr]}_pt${pt1min[rr]}-${pt1max[rr]}_y${dy[rr]/'.'/'p'}_${label}
    [[ $RUN_SAVETREE -eq 1 ]] && { build/ddbar_savetree "$inputdata" $inputmc "$inputmixdata" $outputdir ${pt1min[rr]} ${pt1max[rr]} ${pt2min[rr]} ${pt2max[rr]} ${dy[rr]} ${cmin[rr]} ${cmax[rr]} $label $eff $swapmc $mctree; }
    [[ $RUN_2DFIT -eq 1 ]] && { build/ddbar_2dfit "rootfiles/$outputdir/$dpairtree" "rootfiles/$outputdir/${swaptree}" $outputdir ${pt1min[rr]} ${pt1max[rr]} ${pt2min[rr]} ${pt2max[rr]} ${dy[rr]} ${cmin[rr]} ${cmax[rr]} $label; }
done
