#!/bin/bash

# This script preprocesses data for the TM2D (2 dimensional Tibaldi Molteni) blocking algorithm.
# v0.3: Add CERA20-c
# v0.4: Add ERA40 (for Mik)
# v0.5: ?
# v0.6: Added more datasets

echo "Starting programm TMblock_v0.6 `date +%d:%T.%3N` PID:$$"

declare remapping=true     # [true/false] Remapping the data, if yes: res must be given --> see cdo for naming conventions
res=r180x91 #T63           # r180x91 stands for 2 degree regular lon-lat file

declare download=true         # [Y/N] Copy data from climstor
declare preproc=true          # [Y/N] Preprocessing
declare tmtrack=true          # [Y/N] Compute blockings
declare pproc=true            # [Y/N] Run post processing


declare datatype="era5"    # Type of dataset (must be identical to folder and filename!
declare startyear=2000
declare endyear=2017

declare variable="z500"    # Variable name as given in 

declare persistence=20     # Tuning parameter for TM2D: Minimum persistence (in timesteps)
declare overlap=0.7        # Tuning parameter for TM2D: Minimum overlap between two timesteps
declare tmy=15                     # Tuning parameter for TM2D: Latitude difference between the northern and central / central and southern latitude

# Path to required programs
declare f90_tm2d=tm2d
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"`nf-config --prefix`/lib"

 #Directories
rootdir=`pwd`
if [ ! -d ${rootdir}/workdir ] ; then mkdir ${rootdir}/workdir ; fi

case ${res} in
   T63)
     spherical="TRUE"
     ;;
   *)
     spherical="FALSE"
     ;;
esac # 




case ${datatype} in
   20CR|20cr|20CRv2|20crv2)
      basedir="/storage/home/rohrer/20cr/${variable}"
      workdir=${rootdir}/workdir/20cr ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="20cr"
      list=( 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041 042 043 044 045 046 047 048 049 050 051 052 053 054 055 056 )
      outdir="${outdir}/20cr"
      ;;
   20CRv2c|20crv2c)
      basedir="/storage/home/rohrer/20crv2c/${var2}"
      workdir=${rootdir}/workdir/20crv2c ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="20crv2c"
      list=( 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041 042 043 044 045 046 047 048 049 050 051 052 053 054 055 056 )
      outdir="${outdir}/20crv2c"
      ;;
   CERA|CERA-20C|CERA20C|CERA-20c|CERA20c|cera20c|cera-20c)
      basedir="/storage/home/rohrer/Data/cera20c/${var2}"
      workdir="${rootdir}/workdir/cera20c" ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="cera20c"
      list=( 2 )
      outdir="${outdir}/cera20c"
      ;;
   ERA20CM|ERA-20CM|era20CM|ERA-20cm|era-20cm|era20cm)
      basedir="/storage/home/rohrer/Data/era20cm/${var2}"
      workdir="${rootdir}/workdir/era20cm" ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="era20cm"
      list=( 0 1 2 3 4 5 6 7 8 9 )
      outdir="${outdir}/era20cm"
      ;;
   ERA20c|era20c|ERA20C|era20C|ERA-20c|ERA-20C|era-20c|era-20C)
      basedir="/storage/home/rohrer/Data/era20c/${var2}"
      workdir=${rootdir}/workdir/era20c ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="era20c"
      list=( 000 )
      outdir="${outdir}/era20c"
      ;;
   ERAint|eraint|ERA-interim|ERAi|erai|ERA-int|ERA-INT)
      basedir="/storage/home/rohrer/Data/eraint/${var2}"
      workdir=${rootdir}/workdir/eraint
      datatype="eraint"
      list=( 000 )
      outdir="${outdir}/eraint"
      ;;
   ERA-40|ERA40|era-40|era40|Era-40|Era40)
      basedir="/storage/home/rohrer/Data/era40/${var2}"
      workdir=${rootdir}/workdir/era40
      datatype="era40"
      filename="era40"
      list=( 000 )
      outdir="${outdir}/era40"
      ;;
   ERA-presat|ERA-PRESAT|erapresat|ERApresat|era-presat)
      basedir="/storage/home/rohrer/Data/erapresat/${var2}"
      workdir=${rootdir}/workdir/erapresat
      datatype="erapresat"
      filename="erapresat"
      list=( 000 )
      outdir="${outdir}/erapresat"
      ;;
   JRA|JRA55|jra|jra55|JRA-55|jra-55)
      basedir="/storage/home/rohrer/Data/jra55/${var2}"
      workdir=${rootdir}/workdir/jra55
      datatype="jra55"
      filename="jra55"
      list=( 000 )
      outdir="${outdir}/jra55"
      ;;
   ccc|CCC|CCC400|ccc400)
      #basedir="/mnt/climstor/ccc400"
      basedir="/storage/home/rohrer/ccc400/rawout/Z500"
      workdir=${rootdir}/workdir/ccc400 ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="ccc"
      #list=( 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 
      list=(016 017 018 019 020 021 022 023 024 025 026 027 028 )
     #list=( 029 030 103 )
      list=( 101 102 )
      outdir="${outdir}/ccc400"
      ;;
   ccc_novolc|cccnovolc|ccc400_novolc|ccc400novolc)
      basedir="/storage/home/rohrer/Data/ccc400_novolc/Z500"
      workdir=${rootdir}/workdir/ccc400_novolc ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="ccc_novolc"
      #list=( 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 ) 
      #list=( 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 )
      list=( 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 )
      ;;
   ccc_climsst|ccc400_climsst|cccsst|ccc400sst)
      basedir="/storage/home/rohrer/Data/ccc400_climsst/Z500"
      workdir=${rootdir}/workdir/ccc400_climsst ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="ccc_climsst"
      list=( 001 002 003 004 005 006 007 008 009 010 )
      ;;
   MERRA2|MERRA-2|merra2|merra-2|Merra2|Merra-2)
      basedir="/storage/home/rohrer/Data/merra2/${var2}"
      workdir=${rootdir}/workdir/merra2 ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="merra2"
      filename="merra2"
      list=( 000 )
      outdir="${outdir}/merra2"
      ;;
   MERRA|MERRA-1|merra|merra-1|Merra|Merra-1|merra1|Merra1|MERRA1)
      basedir="/storage/home/rohrer/Data/merra/${var2}"
      workdir=${rootdir}/workdir/merra ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="merra"
      filename="merra"
      list=( 000 )
      outdir="${outdir}/merra"
      ;;
   CFSR|cfsr|NCEP-CFSR)
      basedir="/storage/home/rohrer/Data/cfsr/${var2}"
      workdir=${rootdir}/workdir/cfsr ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="cfsr"
      filename="cfsr"
      list=( 000 )
      outdir="${outdir}/cfsr"
      ;;
   NCEPNCAR|NCEP|ncepncar|ncep|Ncep|NCEP-NCAR|Ncep-Ncar|ncep-ncar)
      basedir="/storage/home/rohrer/Data/ncep/${var2}"
      workdir=${rootdir}/workdir/ncep ; if [ ! -d ${workdir} ] ; then mkdir ${workdir} ; fi
      datatype="ncep"
      filename="ncep"
      list=( 000 )
      outdir="${outdir}/ncep"
      ;;
   *)
      echo "Datatype not found... Check your input"
      exit 1
      ;;
esac

case ${variable} in
   z200)
      var2="Z200" ; ccclev=20000
      ;;
   z500)
      var2="Z500" ; ccclev=50000
      ;;
esac # variable

#######################################################################################
#--------------------------------------------------------------------------------------
#                                        Run
#--------------------------------------------------------------------------------------
#######################################################################################

case $datatype in
   "20cr")
      ;;
   "20crv2c") 
      ;;
   "era5")
      work_d=${workdir}
      if [ ! -d ${work_d} ] ; then mkdir ${work_d}{,/raw${var2},/${var2},/blockings,/anom} ; fi
      if [ "${download}" == "Y" ] ; then
         echo "Download and check ERA5 : `date +'%d:%H:%M:%S:%3N'`"
         for year in `seq ${startyear} 1 ${endyear}` ; do
            if [ "${year}" == "${endyear}" ] ; then
               mnmax=${endmon}
            else
               mnmax=12
            fi # find last month available of year
            for mn in `seq -f %02g 1 ${mnmax}` ; do
               inputfile=${basedir}/${datatype}.${var2}.${year}${mn}.nc4z
               outputfile=${workdir}/${var2}/${var2}_${year}${mn}.nc
               # Change latitude form -180,180 to 0,360 so grid is ok for our tools 
               cdo -s sellonlatbox,0,359.99,-90,90 -chname,${var2},Z ${inputfile} ${outputfile}
            done # for month do standardize data
            cdo -s -O -f nc4 -z zip_1 mergetime ${workdir}/${var2}/${var2}_${year}??.nc ${workdir}/${var2}/${var2}_${year}.nc
         done # for year
      fi # if download
	  ;;
   "merra"|"merra2"|"ncep"|"ncep2"|"jra55"|"cfsr"|"era40"|"erapresat"|"eraint")
      work_d=${workdir}
      if [ ! -d ${work_d}/raw${var2} ] ; then mkdir ${work_d}{,/raw${var2},/${var2},/blockings,/anom} ; fi
      if [ "${download}" == "Y" ] ; then
         echo "Download and check ${filename} : `date +'%d:%H:%M:%S:%3N'`"
         for year in `seq ${startyear} 1 ${endyear}` ; do
            echo "Preprocess year ${year}: `date +'%d:%H:%M:%S:%3N'`"
            inputfile=${basedir}/${filename}.${var2}.${year}.nc
            outputfile=${workdir}/${var2}/${var2}_${year}.nc
            echo "${f90_check} ${inputfile} ${outputfile} ${year} -1 F ${var2} ${datatype}"
            # Change latitude form -180,180 to 0,360 so grid is ok for our tools
            cdo -L -s sellonlatbox,0,359.99,-90,90 -chname,${var2},Z ${work_d}/raw${var2}/${datatype}.${var2}_${year}.nc ${workdir}/${var2}/${var2}_${year}.nc
         done # for year
      fi # if download
      ;;
   "cera20c"|"era20c")
      if [ "${download}" == "Y" ] ; then
         for ennr in ${list[@]} ; do
            echo "Download and check member ${ennr}: `date +'%d:%H:%M:%S:%3N'`"
            work_d=${workdir}/${datatype}_${ennr}
            if [ ! -d ${work_d} ] ; then mkdir ${work_d}{,/${var2},/raw${var2},/blockings,/anom} ; fi
            if [ "${datatype}" == "cera20c" ] ; then
               varname="Z"
            else
               varname="z"
            fi # if datatype
            for year in `seq ${startyear} 1 ${endyear}` ; do
               echo "Download ennr: ${ennr} and year: ${year} ; `date +'%d:%H:%M:%S:%3N'`"
               inputfile=${basedir}/${datatype}_${ennr}.${var2}.${year}.nc
               outputfile=${work_d}/${var2}/${var2}_${year}.nc
               cdo -L -s sellonlatbox,0,359.99,-90,90 -chname,${var2},${varname} ${inputfile} ${outputfile}
            done # for year
         done # for ennr
      fi # if download
      ;;
esac


for ennr in ${list[@]} ; do

   case ${datatype} in
      20cr)
        echo $workdir
        work_d=${workdir}/20cr_${ennr} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=20cr_${ennr}
        ;;
      20crv2c)
        work_d=${workdir}/20crv2c_${ennr} ; if [ ! -d ${workdir} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=20crv2c_${ennr}
        ;;
      era20c)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=era20c
        ;;
      eraint)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=eraint
        ;;
      era40)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=era40
        ;;
      jra55)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=jra55
        ;;
      merra)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=merra
        ;;
      merra2)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=merra2
        ;;
      cfsr)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=cfsr
        ;;
      ncep)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=ncep #NCEP-NCAR
        ;;
      ccc|ccc_novolc|ccc_climsst)
        work_d=${workdir}/ccc400_${ennr} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=ccc400_${ennr}
        ;;
      cera20c)
        work_d=${workdir}/cera20c_${ennr} ;  if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=cera20c_${ennr}
        ;;
      era20cm)
        work_d=${workdir}/era20cm_${ennr} ;  if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=era20cm_${ennr}
        ;;
      erapresat)
        work_d=${workdir} ; if [ ! -d ${work_d} ] ; then echo "Where is your input data then?" ; exit 1 ; fi
        typ=erapresat
        ;;
   esac
   enr=$(echo $ennr | sed 's/^0*//') # Remove trailing zeroes

   if [ "${preproc}" = "Y" ] ; then 
      for year in `seq ${startyear} 1 ${endyear}`; do
         if [ "${remapping}" = "Y" ] ; then
            echo "Remap year ${year}: `date +%d:%T.%3N`"
            if [[ "${spherical}" == "TRUE" ]] ; then
               cdo  -b F32 genbil,n180 ${work_d}/${var2}/${var2}_${year}.nc remapweights.nc
               cdo  -b F32 remap,n180,remapweights.nc ${work_d}/${var2}/${var2}_${year}.nc ${work_d}/${var2}/temp.gauss.nc
               cdo  -b F32 sp2gp -sp2sp,63 -gp2sp ${work_d}/${var2}/temp.gauss.nc ${work_d}/${var2}/temp.${var2}_${year}.nc
            else 
               cdo -s remapbil,${res} ${work_d}/${var2}/${var2}_${year}.nc ${work_d}/${var2}/temp.${var2}_${year}.nc
            fi
         fi # if preproc
      done # for year
      echo "Merge years `date +%d:%T.%3N`"
      if [ "${remapping}" = "Y" ] ; then
         cdo -s -a mergetime ${work_d}/${var2}/temp.${var2}_????.nc ${work_d}/${var2}/${var2}.nc
      else 
         cdo -s -a mergetime ${work_d}/${var2}/${var2}_????.nc ${work_d}/${var2}/${var2}.nc
      fi # if remap
      echo "Merge years done `date +%d:%T.%3N`"

   fi # if preproc

   if [ "${tmtrack}" = "Y" ] ; then
      echo "Track blockings start: `date +%d:%T.%3N`"
      if [ "${datatype}" == "ccc_novolc" ] ; then
          dataname=${work_d}/blockings/blocks_${datatype}_${startyear}-${endyear}.nc 
      else 
          dataname=${work_d}/blockings/blocks_${datatype}.nc
      fi
	  # Actual calculation of TM2D
      ${f90_tm2d} --infile=${work_d}/${var2}/${var2}.nc --invar=Z --outfile=${dataname} --mode=TM2D --persistence=${persistence} --overlap=${overlap} --dlat=${tmy}

      if [ "${datatype}" == "ccc_novolc" ] ; then
         mv blockstat.txt ${work_d}/blockings/blockstat_${startyear}-${endyear}.txt
         mv blocktracks.txt ${work_d}/blockings/blocktracks_${startyear}-${endyear}.txt
         mv blockstat_all.txt ${work_d}/blockings/blockstat_all_${startyear}-${endyear}.txt
         mv blocktracks_all.txt ${work_d}/blockings/blocktracks_all_${startyear}-${endyear}.txt
      else 
         mv blockstat.txt ${work_d}/blockings/blockstat.txt
         mv blocktracks.txt ${work_d}/blockings/blocktracks.txt
         mv blockstat_all.txt ${work_d}/blockings/blockstat_all.txt
         mv blocktracks_all.txt ${work_d}/blockings/blocktracks_all.txt
      fi # if ccc_novolc
   fi # if tmtrack
   
   if [ "${pproc}" == "Y" ] ; then
      echo "Postprocessing: `date +%d:%T.%3N`"
      cdo -s -b F32 seasmean ${dataname} ${work_d}/blockings/temp.seasons.nc
      for seas in DJF MAM JJA SON ; do 
         cdo -s selseas,${seas} ${work_d}/blockings/temp.seasons.nc ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc
         if [ "${seas}" == "DJF" ] ; then
            nsteps=`cdo -s ntime ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc`
            endstep=$(( $nsteps - 1 ))
            cdo -s seltimestep,2/${endstep} ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc ${work_d}/blockings/test.nc
            mv ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc ${work_d}/blockings/old.nc
            mv ${work_d}/blockings/test.nc ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc
         fi # if
         cdo -s timmean ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}_${seas}.nc ${work_d}/blockings/blockfreq_${typ}_${seas}.nc
      done
      cdo -s -b F32 yearmean ${dataname} ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}.nc
      cdo -s -b F32 timmean ${work_d}/blockings/blockfreq_${typ}_${startyear}-${endyear}.nc ${work_d}/blockings/blockfreq_${typ}.nc
      rm ${work_d}/blockings/temp.seasons.nc
      rm ${work_d}/blockings/old.nc

   fi # if pproc

done # for ennr

