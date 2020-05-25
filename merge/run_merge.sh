#!/bin/bash

skim=1
ifs=(0) # input directory
jns=(0) # tree

## ifs
inputdirs=(
    /export/d00/scratch/jwang/DntupleRun2018/crab_Dfinder_20200219_HIZeroBias2_Run2017G-17Nov2017_tkpt0p5dls2p5_v2/ # 0
)

########################################
## >>> do not change lines below >>>  ##
########################################

########################################
nts=( "ntDkpi" "empty" "root")
#jns  0      1       2
########################################

tmp=$(date +%y%m%d%H%M%S)
cp merge.C merge_${tmp}.C
g++ merge_${tmp}.C $(root-config --libs --cflags) -g -o merge_${tmp}.exe || exit 1

for i in ${ifs[@]}
do
    [[ $i -lt ${#inputdirs[@]} ]] || break
    inputdir=${inputdirs[i]}
    IFS='/'; subdir=($inputdir); unset IFS;
    request=${subdir[${#subdir[@]}-1]}
    primedir=${inputdir%%${request}*}

    [[ ! -d $inputdir ]] && continue

    ## ======================================== #

    filelist=filelist_${request}.txt
    [[ -f $filelist ]] && {
        # echo "error: filelist $filelist exits. "
        # echo "remove filelist? (y/n):"
        # rewrite=
        # while [[ $rewrite != 'y' && $rewrite != 'n' ]]
        # do
        #     read rewrite
        #     if [[ $rewrite == 'y' ]] ; then { rm $filelist ; } ;
        #     elif [[ $rewrite == 'n' ]] ; then { echo "please change output file name" ; rm merge_${tmp}.exe ; continue ; } ;
        #     else { echo "please input y/n" ; } ; fi ;
        # done
        rm $filelist
    } 

    ls $inputdir/*.root -d > $filelist

    for j in ${jns[@]}
    do
        ntname=${nts[j]}
        set -x
        output=${primedir}/${request}_skimhltBsize_${ntname}.root
        set +x
        willrun=1
        [[ -f $output ]] && {
            echo "error: output $output exits. "
            echo "remove output? (y/n):"
            rewrite=
            while [[ $rewrite != 'y' && $rewrite != 'n' ]]
            do
                read rewrite
                if [[ $rewrite == 'y' ]] ; then { echo "$output removed" ; rm $output ; } ;
                elif [[ $rewrite == 'n' ]] ; then { echo "please change output file name" ; willrun=0 ; } ;
                else { echo "please input y/n" ; } ; fi ;
            done
        }

        [[ $willrun -eq 0 ]] && continue
        [[ ${1:-0} -eq 1 ]] && { ./merge_${tmp}.exe $output $filelist $ntname $skim ; }
    done
done

rm merge_${tmp}.exe
rm merge_${tmp}.C














##

    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1_1031_NoJSON/ # 0
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v2_1031_NoJSON/ # 1
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v2_1031_NoJSON_Run327527_327564/ # 2
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v2_1031_NoJSON_ToComplete/ # 3
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1_1031_NoJSON_ToComplete/ # 4
