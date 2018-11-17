#!/bin/bash

# Sample script to process vanilla 20CRv3 files


declare preprocess=false
declare merge=false
declare tmtrack=false
declare pproc=true

declare remapping=true
declare res=r180x91
declare syear=1836
declare eyear=2015
declare typ="20crv3_ens"

declare f90_tm2d=`pwd`/tm2d
declare pers=20
declare ovlp=0.7

declare rootdir=/oeschgerstor/climatology/reanalyses/20crv3/Z500
declare workdir=~/tank/blocks/workdir ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi

if [ "${remapping}" = "true" ] ; then
   declare work_d=${workdir}/20crv3_ens.${res} ; if [ ! -d ${work_d} ] ; then mkdir ${work_d}{,/Z500,/blocks} ; fi
else
   declare work_d=${workdir}/20crv3_ens ; if [ ! -d ${work_d} ] ; then mkdir ${work_d}{,/Z500,/blocks} ; fi
fi # Include resolution in workdir

if [ "$preprocess" = "true" ] ; then
    echo "Preprocessing start: `date +%d:%T.%3N`"
    for year in `seq ${syear} ${eyear}` ; do
        echo "$year"
        if [ ${year} -le 1980 ] ; then
           vers=v451
        else
           vers=v452
        fi # determine wheter v451 or v452 is necessary

        if [ "${remapping}" = "true" ] ; then
            cdo -L -s -b F32 -remapbil,${res} -sellonlatbox,0,359.99,-90,90 -selhour,0,6,12,18 ${rootdir}/HGT500.${year}_${vers}.nc ${work_d}/Z500/temp.nc
        else 
            cdo -L -s -b F32 -sellonlatbox,0,359.99,-90,90 -selhour,0,6,12,18 ${rootdir}/HGT500.${year}_${vers}.nc ${work_d}/Z500/temp.nc
        fi # if remapping
        ncwa -O -a plev ${work_d}/Z500/temp.nc ${work_d}/Z500/20crv3_ens.${year}.nc
    done # for year
fi # if preprocess

if [ ${merge} = "true" ] ; then
    echo "Merge years start: `date +%d:%T.%3N`"
    cdo -a -L -O mergetime ${work_d}/Z500/20crv3_ens.????.nc ${work_d}/Z500/20crv3_ens.nc
fi # if merge

if [ ${tmtrack} = "true" ] ; then
    echo "Track blockings start: `date +%d:%T.%3N`"
    infile=${work_d}/Z500/20crv3_ens.nc
    outfile=${work_d}/blocks/blocks_20crv3_ens.nc
    # Actual calculation of TM2D
    ${f90_tm2d} --infile=${infile} --invar=HGT --outfile=${outfile} --mode=TM2D --persistence=${pers} --overlap=${ovlp}

    # Move tracks and stats
    mv blockstat.txt ${work_d}/blocks/blockstat.txt
    mv blocktracks.txt ${work_d}/blocks/blocktracks.txt
    mv blockstat_all.txt ${work_d}/blocks/blockstat_all.txt
    mv blocktracks_all.txt ${work_d}/blocks/blocktracks_all.txt
fi # tmtrack

if [ "${pproc}" = "true" ] ; then
    echo "Postprocessing: `date +%d:%T.%3N`"
    cdo -s -b F32 seasmean ${work_d}/blocks/blocks_20crv3_ens.nc ${work_d}/blocks/temp.seasons.nc
    for seas in DJF MAM JJA SON ; do
       cdo -s selseas,${seas} ${work_d}/blocks/temp.seasons.nc ${work_d}/blocks/blockfreq_${typ}_${syear}-${eyear}_${seas}.nc
       if [ "${seas}" == "DJF" ] ; then
            nsteps=`cdo -s ntime ${work_d}/blocks/blockfreq_${typ}_${syear}-${eyear}_${seas}.nc`
            endstep=$(( $nsteps - 1 ))
            # The first and the last year in the seasmean for DJF are not complete (only Dec for last year, only JF for first year --> Remove)
            cdo -s seltimestep,2/${endstep} ${work_d}/blocks/blockfreq_${typ}_${syear}-${eyear}_${seas}.nc ${work_d}/blocks/test.nc
            mv ${work_d}/blocks/blockfreq_${typ}_${syear}-${eyear}_${seas}.nc ${work_d}/blocks/old.nc
            mv ${work_d}/blocks/test.nc ${work_d}/blocks/blockfreq_${typ}_${syear}-${eyear}_${seas}.nc
        fi # if
        cdo -s timmean ${work_d}/blocks/blockfreq_${typ}_${syear}-${eyear}_${seas}.nc ${work_d}/blocks/blockfreq_${typ}_${seas}.nc
    done
    cdo -s -b F32 yearmean ${work_d}/blocks/blocks_20crv3_ens.nc ${work_d}/blocks/blockfreq_${typ}_${syear}-${eyear}.nc
    cdo -s -b F32 timmean ${work_d}/blocks/blockfreq_${typ}_${syear}-${eyear}.nc ${work_d}/blocks/blockfreq_${typ}.nc
    rm ${work_d}/blocks/temp.seasons.nc
    rm ${work_d}/blocks/old.nc
fi # pproc

