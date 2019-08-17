#!/bin/bash

# This script preprocesses data for the TM2D (2 dimensional Tibaldi Molteni) blocking algorithm.
# v0.3: Add CERA20-c
# v0.4: Add ERA40 (for Mik)
# v0.5: ?
# v0.6: Added more datasets
# v0.7: Add ERA5
# v0.8: Add NCEP-DOE
# v0.8: Complete rewrite to PV (only ERA-5 and ERA-interim, 

echo "Starting programm PVblock_v0.8 `date +%d:%T.%3N` PID:$$"

declare remapping=N   # [Y/N] Remapping the data, if yes: res must be given --> see cdo for naming conventions
res=r360x181 #360x181 #T63                # r180x91 stands for 2 degree regular lon-lat file

declare download=Y         # [Y/N] Copy data from climstor
declare splitter=Y         # [Y/N] Split individual members of 20CR
declare preproc=Y          # [Y/N] Preprocessing
declare tmtrack=Y          # [Y/N] Compute blockings
declare pproc=Y            # [Y/N] Run post processing

bltype="TM2D"
datatype="eraint"    
startyear=1979
endyear=2018

persistence=20             # Tuning parameter for TM2D: Minimum persistence of blocking
overlap=0.7                # Tuning parameter for TM2D: Minimum overlap between timesteps

dlat=15                    # Tuning parameter for TM2D: Latitude difference between the northern and central / central and southern latitude
dzp=-10
anom=200
vapv=-1.3
# Programs 
f90_tm2d=/mnt/data/workdir/script/TM2D/tm2d

# In case that Netcdf library is not found
#export LD_LIBRARY_PATH=/usr/lib


 #Directories

rootdir=/mnt/data/workdir/script/TM2D
if [ ! -d ${rootdir} ] ; then echo "Rootdir not existing, where do you run your script anyway??" ; exit 1 ; fi
if [ ! -d ${rootdir}/workdir ] ; then mkdir ${rootdir}/workdir2 ; fi
outdir=/mnt/data/workdir
#if [ ! -d ${outdir} ] ; then mkdir ${outdir} ; fi
case ${res} in
   T63)
     spherical="TRUE"
     ;;
   *)
     spherical="FALSE"
     ;;
esac # 

if [ $remapping == "N" ] ; then
   resfix="orig"
else 
   resfix=${res}
fi # if remapping


case ${bltype} in
   VAPV)
      suffix="VAPV"
      var="VAPV"
      options="--invar=PV --mode=VAPV --vapvmin=${vapv}"
      declare runmean="Y"
      declare anom="Y"
      declare preprocanom="N"
      declare era5anom="Y"
      ;;
   TM2D)
      suffix="TM2D"
      var="Z500"
      options="--invar=Z500 --mode=TM2D --dlat=${dlat} --dzp=${dzp}"
      declare runmean="N"
      ;;
   Z500anom)
      suffix="Z500anom"
      var="Z500"
      options="--invar=Z500 --mode=Z500anom --anom=${anom}"
      declare runmean="Y"
      declare anom="Y"
      declare preprocanom="N"
      declare era5anom="Y"
      ;;
esac # bltype

case ${datatype} in
   ERAint|eraint|ERA-interim|ERAi|erai|ERA-int|ERA-INT)
      basedir="/mnt/data/workdir/eraint/${var}"
      workdir=${rootdir}/workdir/${suffix}.eraint.${resfix}
      datatype="eraint"
      list=( 000 )
      outdir="${outdir}/eraint"
      ;;
   ERA-5|ERA-5|era5|ERA5)
      basedir="/mnt/data/workdir/era5/${var}"
      workdir=${rootdir}/workdir/${suffix}.era5.${resfix}
      datatype="era5"
      filename="era5"
      list=( 000 )
      outdir="${outdir}/era5"
      ;;
   *)
      echo "Datatype not found... Check your input"
      exit 1
      ;;
esac

#######################################################################################
#--------------------------------------------------------------------------------------
#                                        Run
#--------------------------------------------------------------------------------------
#######################################################################################
# Download data
if [ "${datatype}" == "eraint" ] ; then
   work_d=${workdir}
   if [ ! -d ${work_d} ] ; then mkdir ${work_d}{,/rawPV,/PV,/blockings,/anom} ; fi
   if [ "${download}" == "Y" ] ; then
      echo "Download and check ERA-interim : `date +'%d:%H:%M:%S:%3N'`"

      for year in `seq ${startyear} 1 ${endyear}` ; do
         inputfile=${basedir}/eraint.${var}.${year}.nc
         outputfile=${work_d}/rawPV/ERAi_PV_${year}.nc
         #${f90_check} ${inputfile} ${outputfile} ${year} -1 T PV eraint

         # Change latitude form -180,180 to 0,360 so grid is ok for our tools 
         cdo -s sellonlatbox,0,359.99,-90,90 $inputfile ${workdir}/PV/PV_${year}.nc
      done # for year
   fi # if download ERAint
elif [ "${datatype}" == "merra" ] || [ "${datatype}" == "merra2" ] || [ "${datatype}" == "ncep" ] || [ "${datatype}" == "ncep2" ] || [ "${datatype}" == "jra55" ] || [ "${datatype}" == "cfsr" ]|| [ "${datatype}" == "era40" ] || [ "${datatype}" == "erapresat" ] || [ "${datatype}" == "era5" ] ; then
   work_d=${workdir}
   if [ ! -d ${work_d}/rawPV ] ; then mkdir ${work_d}{,/rawPV,/PV,/blockings,/anom} ; fi

   if [ "${download}" == "Y" ] ; then
     echo "Download and check ${filename} : `date +'%d:%H:%M:%S:%3N'`"

     for year in `seq ${startyear} 1 ${endyear}` ; do
         echo "Preprocess year ${year}: `date +'%d:%H:%M:%S:%3N'`"
         inputfile=${basedir}/${filename}.${var}.${year}.nc
         outputfile=${work_d}/rawPV/${datatype}.PV_${year}.nc
         echo "$inputfile PV"
         #${f90_check} ${inputfile} ${outputfile} ${year} -1 T PV ${datatype}
         # Change latitude form -180,180 to 0,360 so grid is ok for our tools 
         cdo -a -L -s sellonlatbox,0,359.99,-90,90 ${inputfile} ${workdir}/PV/PV_${year}.nc
      done # for year
   fi # if download
fi # if datatype


# Download finished --> Start preproc
for ennr in ${list[@]} ; do
   case ${datatype} in
      eraint)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=eraint
        ;;
      era5)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=era5
        ;;
      merra2)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=merra2
        ;;
      cfsr)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=cfsr
        ;;
   esac
   enr=$(echo $ennr | sed 's/^0*//') # Remove trailing zeroes

   if [ "${preproc}" == "Y" ] ; then
      for year in `seq ${startyear} 1 ${endyear}`; do
         if [ "${remapping}" = "Y" ] ; then
            echo "Remap year ${year}: `date +%d:%T.%3N`"
            if [[ "${spherical}" == "TRUE" ]] ; then
               cdo  -b F32 genbil,n180 ${work_d}/PV/PV_${year}.nc remapweights.nc
               cdo  -b F32 remap,n180,remapweights.nc ${work_d}/PV/PV_${year}.nc ${work_d}/PV/temp.gauss.nc
               cdo  -b F32 sp2gp -sp2sp,63 -gp2sp ${work_d}/PV/temp.gauss.nc ${work_d}/PV/temp.PV_${year}.nc
            else
               cdo -s remapbil,${res} ${work_d}/PV/PV_${year}.nc ${work_d}/PV/temp.PV_${year}.nc
            fi
         fi # if preproc
      done # for year
      echo "Merge years `date +%d:%T.%3N`"
      if [ "${remapping}" = "Y" ] ; then
         cdo -O -s -a mergetime ${work_d}/PV/temp.PV_????.nc ${work_d}/PV/PV.nc
      else
         cdo -O -s -a mergetime ${work_d}/PV/PV_????.nc ${work_d}/PV/PV.nc
      fi # if remap
      echo "Merge years done `date +%d:%T.%3N`"

   fi # if preproc


   if [ "${runmean}" == "Y" ] && [ "${preproc}" == "Y" ] ; then
      echo "Running mean `date +%d:%T.%3N`"
      # #define leapyear function
      leapyear () {
         if [ $[$1 % 4] -eq 0 ] && [ $[$1 % 100] -ne 0 ] || [ $[$1 % 400] -eq 0 ]; then
            tpy=1464
            leapyear="yes"
         else
            tpy=1460
            leapyear="no"
         fi;
      } # end function

      if [ "${anom}" == "Y" ] ; then
        if [ "${preprocanom}" == "N" ] ; then
            cdo -s runmean,121 ${work_d}/PV/PV.nc ${work_d}/PV/PV_rm121.nc
            if [ "${era5anom}" == "Y" ] ; then
               cdo -s selyear,1981/2010 ${work_d}/PV/PV_rm121.nc ${work_d}/PV/PV_rm121_crop.nc
               cdo -s selyear,1600/2200 -yhourmean -delete,day=29,month=2 ${work_d}/PV/PV_rm121_crop.nc ${work_d}/PV/clim.nc
            else
               cdo -s selyear,1600/2200 -yhourmean -delete,day=29,month=2 ${work_d}/PV/PV_rm121.nc ${work_d}/PV/clim.nc
            fi # era5anom
            cdo -s select,day=28,month=2 ${work_d}/PV/clim.nc ${work_d}/PV/clim28feb.nc
            cdo -s cat -seltimestep,1/236 ${work_d}/PV/clim.nc ${work_d}/PV/clim28feb.nc -seltimestep,237/1460 ${work_d}/PV/clim.nc ${work_d}/PV/clim366.nc
         else
            cp ~/tank/PVblock/VAPV_clim2_365.nc ${work_d}/PV/clim.nc
            cp ~/tank/PVblock/VAPV_clim2_366.nc ${work_d}/PV/clim366.nc
         fi
         for year in `seq ${startyear} 1 ${endyear}` ; do
            echo "Calculate anomaly ${year} : `date +%d:%T.%3N`"
            leapyear ${year}
            if [ leapyear = "yes" ];then
               cdo -s sub -selyear,${year} ${work_d}/PV/PV.nc ${work_d}/PV/clim366.nc ${work_d}/PV/PV_anom_${year}.nc
            else
               cdo -s sub -selyear,${year} ${work_d}/PV/PV.nc ${work_d}/PV/clim.nc ${work_d}/PV/PV_anom_${year}.nc
            fi # if leapyear
         done # for year
         cdo -O -s mergetime ${work_d}/PV/PV_anom_????.nc ${work_d}/PV/PV_anom.nc
         cdo -s runmean,9 ${work_d}/PV/PV_anom.nc ${work_d}/PV/PV_rm.nc
      else
         # Two day running mean  
         cdo -s runmean,9 ${work_d}/PV/PV.nc ${work_d}/PV/PV_rm.nc
      fi # if anom
   else
      ln -s ${work_d}/PV/PV.nc ${work_d}/PV/PV_rm.nc
   fi # if runmean

  dataname=${work_d}/blockings/blocks_${datatype}.nc
  if [ "${tmtrack}" == "Y" ] ; then
      echo "Track blockings start: `date +%d:%T.%3N`"
      echo " ${f90_tm2d} --infile=${work_d}/PV/PV_rm.nc --outfile=${dataname} ${options} --persistence=${persistence} --overlap=${overlap}"
      ${f90_tm2d} --infile=${work_d}/PV/PV_rm.nc --outfile=${dataname} ${options} --persistence=${persistence} --overlap=${overlap}
    
      mv blockstat.txt ${work_d}/blockings/blockstat.txt
      mv blocktracks.txt ${work_d}/blockings/blocktracks.txt
      mv blockstat_all.txt ${work_d}/blockings/blockstat_all.txt
      mv blocktracks_all.txt ${work_d}/blockings/blocktracks_all.txt
   fi # if tmtrack

   if [ "${pproc}" == "Y" ] ; then
      echo "Postprocessing: `date +%d:%T.%3N`"
      cdo -f nc -s -b F32 seasmean -gtc,0.5 ${dataname} ${work_d}/blockings/temp.seasons.nc
      for seas in DJF MAM JJA SON ; do
         cdo -f nc -s selseas,${seas} ${work_d}/blockings/temp.seasons.nc ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc
         if [ "${seas}" == "DJF" ] ; then
            nsteps=`cdo -s ntime ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc`
            endstep=$(( $nsteps - 1 ))
            cdo -s seltimestep,2/${endstep} ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc ${work_d}/blockings/test.nc
            mv ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc ${work_d}/blockings/old.nc
            mv ${work_d}/blockings/test.nc ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc
         fi # if
         cdo -s timmean -gtc,0.5 ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc ${work_d}/blockings/blockfreq_${typ}_${seas}.nc
      done
      cdo -f nc -s -b F32 yearmean -gtc,0.5 ${dataname} ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}.nc
      cdo -s -b F32 timmean ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}.nc ${work_d}/blockings/blockfreq_${typ}.nc
      rm ${work_d}/blockings/temp.seasons.nc
      rm ${work_d}/blockings/old.nc

   fi # if pproc

done # for ennr

