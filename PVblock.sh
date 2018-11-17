#!/bin/bash

# Compute the blocks based on the 500-150 hPa vertically averaged potential vorticity (VAPV). The algorithm is very similar to Schwierz et al. (2004). 
# The script contains several parts:
# Definition part: Describes which datasets, define the threshold for the blocking algorithm
# Preprocessing: 
# Compute anomalies
# Compute blocks
# Postprocessing


echo "Starting programm PVblock_v0.8 `date +%d:%T.%3N` PID:$$"


dlist=( f.e122.F_1850-PDAY_CAM5.f09_g16.climatology.oFaIsF-io192.000
f.e122.F_1850-PDAY_CAM5.f09_g16.climatology.oFaIsF-io192.001
f.e122.F_1850-PDAY_CAM5.f09_g16.climatology.oFaIsF-io192.002 )



declare syear=2009 # First year of dataset (2006 or 2009 for experiments, 1979 for climatologies
declare eyear=2016 # Last year of dataset (2016 for experiments, 2005 for climatologies)

declare preprocess=true
declare runmean=true
declare anom=true
declare precomputedanom=true    # If false, the anomaly is calculated for thsi dataset
declare anomtype="oFaIsF"       # If precomputedanom is taken, indicate here which one
declare anomsyea=1982 # Start year of ANOMALY (do not mix of syear!)
declare anomeyea=2005 # End year of ANOMALY 
declare tmtrack=true
declare pproc=true

# You may change the thresholds that define a block, normally just leave as is

persistence=20             # Tuning parameter for TM2D: Minimum persistence of blocking
overlap=0.7                # Tuning parameter for TM2D: Minimum overlap between timesteps
vapvmin=-1.3               # Tuning parameter for TM2D: PV minimum threshold

# Programs 
tm2d=/home/marco/TM2D/bin/tm2d

 #Directories
rootdir=`pwd`
if [ ! -d ${rootdir} ] ; then echo "Rootdir not existing, where do you run your script anyway??" ; exit 1 ; fi
if [ ! -d ${rootdir}/workdir ] ; then mkdir ${rootdir}/workdir ; fi
indir=/landclim1/wehrlik/ExtremeX/${dset}/atm/hist



# Start with calculations
for dset in ${dlist[@]} ; do
   echo "Start ${dset}: `date +%d:%T.%3N`"
   if [ ! -d ${rootdir}/workdir/${dset} ] ; then mkdir ${rootdir}/workdir/${dset} ; fi
   if ( ${precomputedanom} ) then
      workdir=${rootdir}/workdir/${dset}/${anomtype}.${anomsyea}-${anomeyea}
      if [ ! -d ${rootdir}/workdir/${dset}/${anomtype}.${anomsyea}-${anomeyea} ] ; then mkdir ${rootdir}/workdir/${dset}/${anomtype}.${anomsyea}-${anomeyea}{,/PV,/indat,/blocks} ; fi
   else
      workdir=${rootdir}/workdir/${dset}
      if [ ! -d ${rootdir}/workdir/${dset}/ownanom ] ; then mkdir ${rootdir}/workdir/${dset}/ownanom{,/PV,/indat,/blocks} ; fi
   fi # if takeeraintanom

### Get data - extract PV field form model output and make yearly files
   if ( ${preprocess} ) then
      echo "Preprocessing of ${dset}: `date +%d:%T.%3N`"
      for year in `seq ${syear} 1 ${eyear}` ; do
         n=0
         filetype=`echo ${dset:(-7):3}`
         if [ ${filetype} -eq 192 ] ; then
            for file in `ls ${indir}/*.cam.h2.${year}-*0_pp.nc` ; do
               n=$((n+=1))
               # May be insert option to remap, but I don't think we need it
               cdo -s selvar,APV ${file} ${workdir}/indat/PV${year}.${n}.nc
            done # for file
            cdo -s -a -O mergetime ${workdir}/indat/PV${year}.*.nc ${workdir}/indat/PV${year}.nc
         else
            cdo -s selvar,APV ${indir}/*.cam.h2.${year}-*0_pp.nc ${workdir}/indat/PV${year}.nc
         fi # 
      done # for year    
      cdo -s -a -O mergetime ${workdir}/indat/PV????.nc ${workdir}/indat/${dset}.PV.nc
   fi # if preprocess

### Compute anomalies (either take preexisting climatology or derive climatology from data.  - 
   if ( ${runmean} ) then
      echo "Running mean ${dset} `date +%d:%T.%3N`"
      if ( ${anom} ) then
         if ( ${precomputedanom} ) then
            nn=`ls ${rootdir}/climatologies/${anomtype}/PVclim.${anomtype}.${anomsyea}-${anomeyea}.nc | wc -l`
            if [ ${nn} -eq 1 ] ; then
               cp ${rootdir}/climatologies/${anomtype}/PVclim.${anomtype}.${anomsyea}-${anomeyea}.nc ${workdir}/PV/clim.nc
               touch ${workdir}/PV/PVclim.${anomtype}
            else
               echo "ERROR: Climatology not found ${anomtype}"
               stop
            fi # if nn
         else
            cdo -s runmean,121 -selyear,${syear}/${eyear} ${workdir}/indat/${dset}.PV.nc ${workdir}/PV/PV_rm121.nc
            cdo -s yhourmean -delete,day=29,month=2 ${workdir}/PV/PV_rm121.nc ${workdir}/PV/clim.nc
         fi # if takeeraintanom
         for year in `seq ${syear} 1 ${eyear}` ; do
            echo "Calculate anomaly ${year} : `date +%d:%T.%3N`"
            cdo -s sub -selyear,${year} ${workdir}/indat/${dset}.PV.nc ${workdir}/PV/clim.nc ${workdir}/PV/PV_anom_${year}.nc
         done # for year
         cdo -O -s mergetime ${workdir}/PV/PV_anom_????.nc ${workdir}/PV/PV_anom.nc
         cdo -s runmean,9 ${workdir}/PV/PV_anom.nc ${workdir}/PV/PV_rm.nc
      else
         # Two day running mean  
         cdo -s runmean,9 ${workdir}/indat/${dset}.PV.nc ${workdir}/PV/PV_rm.nc
      fi # if anom
   fi # if runmean

### Track blockings
   dataname=${workdir}/blocks/blocks_${dset}.nc
   if ( ${tmtrack} ) then
      echo "Track blockings start ${dset}: `date +%d:%T.%3N`"
      ${tm2d} --infile=${workdir}/PV/PV_rm.nc --invar=APV --outfile=${dataname} --mode=VAPV --persistence=${persistence} --overlap=${overlap} --vapvmin=${vapvmin}
      mv blockstat.txt ${workdir}/blocks/blockstat.txt
      mv blocktracks.txt ${workdir}/blocks/blocktracks.txt
      mv blockstat_all.txt ${workdir}/blocks/blockstat_all.txt
      mv blocktracks_all.txt ${workdir}/blocks/blocktracks_all.txt
   fi # if ttruek

### Postprocessing --> Create yearly and seasonal mean files
   if ( ${pproc} ) then
      typ=${dset}
      echo "Postprocessing start ${dset}: `date +%d:%T.%3N`"
      cdo -s -b F32 seasmean -gtc,0.5 ${dataname} ${workdir}/blocks/temp.seasons.nc
      for seas in DJF MAM JJA SON ; do
         seasfile="${workdir}/blocks/blockfreq_${typ}_${syear}-${eyear}_${seas}.nc"
         cdo -s selseas,${seas} ${workdir}/blocks/temp.seasons.nc ${seasfile}
         if [ "${seas}" == "DJF" ] ; then
            nsteps=`cdo -s ntime ${seasfile}`
            endstep=$(( $nsteps - 1 ))
            cdo -s seltimestep,2/${endstep} ${seasfile} ${workdir}/blocks/test.nc
            rm ${seasfile}
            mv ${workdir}/blocks/test.nc ${seasfile}
         fi # if seas
         cdo -s timmean ${seasfile} ${workdir}/blocks/blockfreq_${typ}_${seas}.nc
      done # for seas
      cdo -s -b F32 yearmean -gtc,0.5 ${dataname} ${workdir}/blocks/blockfreq_${typ}_${syear}-${eyear}.nc
      cdo -s -b F32 timmean ${workdir}/blocks/blockfreq_${typ}_${syear}-${eyear}.nc ${workdir}/blocks/blockfreq_${typ}.nc
      rm ${workdir}/blocks/temp.seasons.nc
   fi # if pproc

done # for ennr


