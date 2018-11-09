program tm2d

! This is a 2D blocking index using a similar approach as Scherrer et al. (2006).
! Prerequisites: - NetCDF library (libnetcdff
!                - An absolute timeaxis --> use cdo -a to change that
 
! v0.1     : First definition of GPH
! v0.2     : Implement more TM like definition
! v0.5     : Implemented 2D contour searching and overlapping routine (after two failures)
! v0.6     : Implemnt 3D contour searching (not in a subroutine, we do all contours at once for memory/computational reasons
! v0.7     : Blocktrack/blockstat and temporal persistence
! v0.7.1   : Added stuff, reduce memory footprint, debugging
! v0.7.2   : Added stuff, reduce memory footprint, debugging, performance
! v0.7.3   : I/O optimizations, bugfixing (Max dz was wrong)
! v0.8     : Insert gettime, update outputfiles accordingly, use of CDI for writing outputfiel
! v0.8.1   : Revert back to netcdf library (CDI is not so well documented and works not like intended...
!            Start rewriting to make use of time slices --> We only allocate an array for one year and constantly refill this domain until finished
! v0.8.2   : Rewritte/Debug everything upto Contour3D
! v0.8.4   : Rewrite contour3D and debug everything
! v0.8.5   : Code clean up and bug fixes
! v0.8.6   : Code clean up (time steps) and bug fixes (time axis and blocktracks latitude was reversed)
! v0.9.0   : Rewrite for more convenient command line interface, added Netcdf4 compression, now incorporates TM2D, Z500' and VAPV

! Marco Rohrer, Oct 2016 (marco.rohrer@giub.unibe.ch)


USE netcdf

IMPLICIT NONE


LOGICAL                                    :: nopersistence=.FALSE. !.TRUE. !FALSE.
LOGICAL                                    :: nooverlap=.FALSE. !.TRUE. !.FALSE.
LOGICAL                                    :: compression=.TRUE.                               ! Zipped NetCDF-4 file or not?
LOGICAL                                    :: labeling=.FALSE.                                ! Switch whether blocked areas are labelled or binary

! Argument I/O
INTEGER                                    :: narg, i                                          ! Number and counter of input arguments
CHARACTER(128)                             :: arg                                              ! Input argument that will be parsed
INTEGER                                    :: errors=0                                         ! Number of errors

! NetCDF I/O
CHARACTER(1000)                            :: cfoutcomment, cftimeunit, cfcalendar             ! Variables to store NetCDF comments and time units
CHARACTER(200)                             :: infile="unknown", invarname="unknown"            ! Name of input file and input variable (mandatory!)
CHARACTER(200)                             :: inlat, inlon, inver, intime                      ! Names of the incoming file, the variable name that should be read and the dimensions
CHARACTER(200)                             :: outfile="unknown", mode="unknown"                ! Name of outgoing file, calculation mode (mandatory!)
CHARACTER(200)                             :: cflevunit, cflevname                             ! Level attributes
CHARACTER(24),DIMENSION(:),ALLOCATABLE     :: dimnames,varnames                                ! Names of dimensions and variables
INTEGER,DIMENSION(:), ALLOCATABLE          :: lendims                                          ! Dimensions lengths
INTEGER                                    :: ndims,nvars, nglobatt                            ! Number if dimension, variables and global attributes
INTEGER                                    :: ncid, ncido                                      ! NetCDF ID for incoming and outgoing file
INTEGER                                    :: istat                                            ! Status of NetCDF operation
INTEGER                                    :: LatDimID=-1,LonDimID=-1, VerDimID=-1,TimeDimID=-1
INTEGER                                    :: TimeVarID,LonVarID,LatVarID,VerVarID,DatVarID
INTEGER                                    :: LatDimIDo,LonDimIDo,TimeDimIDo
INTEGER                                    :: LatVarIDO,LonVarIDo,TimeVarIDo,DatVarIDo
INTEGER                                    :: nLats,nLong,nTime,nVert
LOGICAL                                    :: exlev=.TRUE.
INTEGER,DIMENSION(:), ALLOCATABLE          :: vertical
REAL,DIMENSION(:), ALLOCATABLE             :: latitude, longitude
REAL(8), DIMENSION(:), ALLOCATABLE         :: times
LOGICAL                                    :: reverselat=.FALSE.


! Arraya
INTEGER                                    :: tlen=720  ! Note, a block cannot be longer than the amount of timesteps declared here!
REAL,DIMENSION(:,:),ALLOCATABLE            :: arr2d     ! Helper array to read/write one time slice
REAL,DIMENSION(:,:,:), ALLOCATABLE         :: arr       ! Array with input data
INTEGER,DIMENSION(:,:,:), ALLOCATABLE      :: arr2,arr3 ! Arrays to calculate the blocking index, two because we have the overlap and the connect3D part which run at the same time and would "destroy" data
REAL,DIMENSION(:,:,:,:), ALLOCATABLE       :: arrin     ! Helper array if the input file has 4 dimension (lon,lat,lev,time)

! Random stuff
CHARACTER(100)                             :: version
CHARACTER(5)                               :: cpers,cover,ctmy,clogf ! Character persistence, overlap and  --> If given as argument
! Count variables
INTEGER                                    :: ii,jj,tt  ! Count variables for lon,lat,time
INTEGER                                    :: ix,iy,it  ! Position of grid point in connect3D/2D, it is the relative time in the array
INTEGER                                    ::       iu  ! iu is the absolute time in the array
INTEGER                                    :: ik,jl,tm  ! Count variables to search for neighbouring points in connect3D/2D
INTEGER                                    :: tx        ! Time index within matrix
INTEGER                                    :: t1,t2,t3  ! The three timesteps for overlap (we need three indices as we can have transistions from the end of the timeaxis to the beginning

! TM2D stuff
REAL                                       :: overlap=0.7     ! TUNING: Indicates the overlap between two time steps
INTEGER                                    :: persistence=20  ! TUNING: Indicates the minimal number of time steps a blocking must endure

REAL                                       :: dx,dy           ! spacing in x (lon), y (lat)
INTEGER                                    :: nstep           ! Vertical steps needed per hemisphere
INTEGER                                    :: minlat=35       ! Minimum equatorward latitude, where algorithm calcs
INTEGER                                    :: tmy=15          ! Tibaldi-Molteni GPH check latitude difference
INTEGER                                    :: dzp=-10         ! Blocking criterion pole   (Zpo-Zmi)/deg
INTEGER                                    :: dze=0           ! Blocking criterion equator(Zmi-Zeq)/deg
REAL                                       :: anom=200.       ! Anomaly criterion (for Z500anom)
REAL                                       :: vapv=-1.3       ! Anomaly criterion (for VAPV)
INTEGER                                    :: tmo             ! Offset resulting from dy and tmy
REAL                                       :: sy,ey           ! Start x, Stop x
INTEGER                                    :: nclu            ! Cluster number
LOGICAL                                    :: nh=.FALSE.      ! Compute northern hemisphere if available
LOGICAL                                    :: sh=.FALSE.      ! Compute southern hemisphere if available
INTEGER                                    :: label2d         ! To discriminate 2D contours
INTEGER                                    :: label=0         ! To discriminate contours (3D)
INTEGER,DIMENSION(:,:), ALLOCATABLE        :: carr            ! Domain to be checked whether it is connected or not
INTEGER,DIMENSION(:,:,:),ALLOCATABLE       :: arrover         ! Contains three timesteps handed over to get overlap
INTEGER,DIMENSION(:),ALLOCATABLE           :: iiv,jjv,ttv,ttva! Vectors which holds the position of all blocked grid points of a block (ttv relative time; ttva absolute time)
REAL,DIMENSION(:),ALLOCATABLE              :: dzv             ! Vector which holds the blocking anomaly per blocked grid point
INTEGER                                    :: ngp,ngpc        ! The number of blocked grid points, counter variable which checks whether blocked neighbours occur for each blocked grid point
REAL                                       :: dzmax           ! Maximum blocking anomaly of whole block

! blockstat/tracks stuff
LOGICAL                                    :: logging=.TRUE.       ! Produce blocktracks and blockstat?
INTEGER                                    :: year,mon,day,hour    ! Date for log file
INTEGER                                    :: dzmaxt,dzmaxx,dzmaxy ! Position of maximum anomaly of block
INTEGER                                    :: nn=0                 ! Helper variable to get position in btracks
REAL,DIMENSION(:,:),ALLOCATABLE            :: bstat,btrack         ! Arrays for blocking stats and tracks
REAL,DIMENSION(:),ALLOCATABLE              :: dzmaxv               ! Vector that holds the maximum blocking anomaly per time step
INTEGER,DIMENSION(:),ALLOCATABLE           :: dzmaxxv,dzmaxyv    
REAL(8),DIMENSION(:),ALLOCATABLE           :: cossum,sinsum        ! Get all sine and cosine values of all blocked gridpoints per timesteps
REAL,DIMENSION(:),ALLOCATABLE              :: latsum
INTEGER,DIMENSION(:),ALLOCATABLE           :: ngpsum
REAL(8)                                    :: pi=4.0*ATAN(1.D0)    ! equals pi=3.14...
!! ====================================================================================

INTERFACE ! Interface for gettime, so we don't need pass nTime for allocation, more educational purpose than any sense in it...
  SUBROUTINE gettime(times,tt,cftimeunit,year,mon,day,hour)
      IMPLICIT NONE

      INTEGER, INTENT(IN)                           :: tt
      REAL(8), DIMENSION(:), INTENT(INOUT)          :: times
      CHARACTER(1000),INTENT(IN)                    :: cftimeunit
      INTEGER,INTENT(OUT)                           :: year,mon,day,hour

   END SUBROUTINE gettime
END INTERFACE


version='v0.9.0'

narg=iargc()
IF (narg.EQ.0) THEN
   PRINT*,"ERROR: Not enough parameters given, at least mode, infile, invar, outfile must be given"
   PRINT*,"Example: ./tm2d --infile=Z500.nc --invar=Z --outfile=blocks.nc --mode==TM2D"
   STOP
ENDIF

PRINT*,"*****************************************************************"

DO i=1,narg
   CALL getarg(i,arg)
   SELECT CASE(arg(1:5)) ! Could be nicer but it works ¯\_(ツ)_/¯
      CASE ('-v   ','--ver')
        PRINT*,"Version 0.9 - Oct 2018, marco.rohrer@giub.unibe.ch"
        STOP
      CASE ('-h   ','--hel')
        PRINT*,"Usage: ./tm2d --infile='Z500.nc' --invar='Z' --outfile='blocks.nc' --mode=='TM2D'"
        PRINT*,"Use --morehelp for a more thourough explanation of the algoritm"
        PRINT*,"Options and explanation:"
        PRINT*,"--morehelp [Print even more help]"
        PRINT*,"--infile=CHAR [Indicate input file name - REQUIRED]"
        PRINT*,"--invar=CHAR [Indicate required variable name of infile - REQUIRED]"
        PRINT*,"--outfile=CHAR [Indicate output file name - REQUIRED]"
        PRINT*,"--mode=CHAR [Indicate blocking mode: {TM2D,Z500anom,VAPV}]"
        PRINT*,"--persistence=INT [# timesteps for persistence, required for TM2D,Z500anom,VAPV][20]"
        PRINT*,"--overlap=FLOAT [required spatial overlap, required for TM2D,Z500anom,VAPV][0.7]"
        PRINT*,"--dlat=INT [degrees between central and N/S grid point for TM2D [15]"
        PRINT*,"--dzp=INT [required gph gradient towards pole [m/°lat]for TM2D [-10]"
        PRINT*,"--dze=INT [required gph gradient towards equator [m/°lat] for TM2D [0]"
        PRINT*,"--minlat=INT [minimum latitude towards equator for TM2D [35]"
        PRINT*,"--anom=INT [minimum geopotential height anomaly [m] for Z500anom [200]"
        PRINT*,"--vapvmin=FLOAT [minimum VAPV treshold [PVU] for VAPV [-1.3]"
        PRINT*,"-nologging [deactivate logging of block characteristics, increases performance]"
        PRINT*,"-nooverlap [deactivate overlap criterium]"
        PRINT*,"-nopersistence [deactivate persistence criterium]"
        PRINT*,"-nc3 [creates an uncompressed NetCDF3 file instead of a zipped NetCDF4 file"
        PRINT*,"-labels [activate ascencding numbering of blocks instead of binary field]"
        STOP
      CASE ('-m   ','--mor')
        WRITE(*,'(A446)') "This program detects and tracks blockings. So far three different block algorithms are implemented that &
             &were used in Rohrer et al. (in prep.). They all operate similarly by first finding areas that fulfil certain requi&
             &rements that qualify them as prototype blocks. Thereafter these 'protoblock areas' are tracked over time and tagge&
             &d as blocks if the protoblock area surpasses an areal overlap test and is long enough (the persistence criterium)."
        WRITE(*,'(A354)') "For all three options the --persistence and the --overlap criterion are valid. Their standard setting&
                & are 0.7 for the overlap and 20 time steps for persistence. For 6-hourly data, 20 time steps correspond to 5 da&
                &ys, that is often used as a threshold for block duration. However these standard settings can be overwritten wi&
                &th --persistence and --overlap."
        WRITE(*,'(A22)') "The three options are:"
        WRITE(*,'(A217,I2,A147,I3,A67,I3,A95)') "a) TM2D:  Geopotential height reversal are detected with the same&
               & approach as Tibaldi and Molteni (1990). See also Tibaldi et al. (1994) for Southern Hemisphere. The&
              & latitudes for the calculation of the gradients are ",tmy,"degrees apart. This can be changed via the&
             & parameter --dlat. The two criteria the detect a block are 1) a meridional gradient that is larger than ",dzp,&
             &" towards the pole and 2) a meridional gradient that is larger than ",dze," towards the equator. These two&
             &  parameter can be changed using the --dzp and --dze parameters."
        WRITE(*,'(A428)') "b) DG83: Positive anomalies in the geopotential height field are detected as protoblocks&
               &  according to the definition in Dole and Gordon (1983). For that matter geopotential height anomalies&
               & should be provided with the --infile option as the algorithm does NOT calculate anomalies on its own.&
               & In its standard setting, geopotential height anomalies must be larger than 200m. This behaviour can be&
               & changed by the --anom parameter."
        WRITE(*,'(A382)') "c) VAPV: Defines blocks as negative anomalies in the vertically averaged potential vorticity&
               & (VAPV) field in the Northern Hemisphere. For easier calculation, the Southern Hemisphere is multiplied by&
               & (*-1) to reverse the effect of the change in sign of the Coriolis parameter => Anomalies in the Southern&
               & Hemisphere must be multiplied by (*-1) for a correct interpretation of values!"
        WRITE(*,*) ""
        WRITE(*,'(A91)') "The program has more options to customize the output [The standard setting is in brackets]:"
        WRITE(*,'(A103)') "INT stands for integer, FLOAT for floating point number, CHAR for a string, BOOL for a boolean&
               & variable"
        WRITE(*,'(A103)') "--logging=BOOL : Use this option to deactivate the log files that saves statistics about the&
               & blocks [T]"
        WRITE(*,'(A105)') "--infile=CHAR: REQUIRED, use this parameter to indicate the file that is read in, e.g.&
               & /home/user/Z500.nc"
        WRITE(*,'(A102)') "--invar=CHAR : REQUIRED, use this parameter to indicate which variable in infile is read on,&
               & e.g. Z500"
        WRITE(*,'(A125)') "--outfile=CHAR : REQUIRED, use this parameter to indicate the file to which the blocks are&
               & written, e.g. /home/user/blocks.nc"
        STOP
   CASE ('-logg')
        logging=.TRUE.
   CASE ('--inf')
        infile=TRIM(arg(10:128))
   CASE ('--inv')
        invarname=TRIM(arg(9:128))
   CASE ('--out')
        outfile=TRIM(arg(11:128))
   CASE ('--mod')
        mode=TRIM(arg(8:32))
   CASE ('--ove')
        READ(arg(11:14),'(F4.0)') overlap
   CASE ('--per')
        READ(arg(15:17),'(I3)') persistence
   CASE ('--dla')
        READ(arg(8:9), '(I2)') tmy
   CASE ('--min')
        READ(arg(10:12),'(I3)') minlat
   CASE ('--dzp')
        READ(arg(7:9), '(I3)') dzp
   CASE ('--dze')
        READ(arg(7:9), '(I3)') dze
   CASE ('--ano')
        READ(arg(8:13), '(I6)') anom
   CASE ('--vap')
        READ(arg(11:16), '(F6.0)') vapv
   CASE ('-nolo')
        logging=.FALSE.
   CASE ('-noov')
        PRINT*, "Attention overlap criterion disabled, overlap set to 0.0"
        overlap=0.0
   CASE ('-nope')
        PRINT*, "Attention persistence criterion disabled, persistence set to 0.0"
        persistence=0
   CASE ('-nc3 ')
        compression=.FALSE.
   CASE ('-labe')
        labeling=.TRUE.
   CASE DEFAULT
        WRITE(*,*) "Attention argument not found:",arg
   END SELECT
ENDDO

! Check if we got the minimum required input
IF (mode=="unknown") THEN
   PRINT*,"--> No mode selected [available: TM2D, VAPV, Z500anom]. Please indicate with --mode=TM2D. Aborting"
   errors=errors+1
ENDIF
IF (infile=="unknown") THEN
   PRINT*,"-->Input file not specified, please indicate with --infile=/path/to/infile. Aborting"
   errors=errors+1
ENDIF
IF (invarname=="unknown") THEN
   PRINT*,"-->Input variable name not specific please indicate with --invar=myvar. &
           &If you don't know, check your input file with ncdump -h. Aborting"
   errors=errors+1
ENDIF
IF (outfile=="unknown") THEN
   PRINT*,"-->Output file not specified, please indicate with --outfile=/path/to/outfile. Aborting"
   errors=errors+1
ENDIF

! Stop, if a mandatory argument is missing
IF (errors.GT.0) STOP


!==================================================================================================================================
!========================== Open NetCDF file 
istat = 0
istat = NF90_OPEN(infile, nf90_nowrite, ncid)
IF (istat /= nf90_NoErr) THEN
   WRITE(*,*) "File "//trim(infile)//" does not exist. Exit"
   STOP
ENDIF

! Get dimensions
istat = 0
istat = NF90_INQUIRE(ncid,ndims,nvars,nglobatt,TimeDimID)
IF (istat/= nf90_NoErr ) WRITE(*,*) "Inquire:", NF90_STRERROR(istat)

ALLOCATE(dimnames(ndims),varnames(nvars),lendims(ndims))
exlev=.FALSE.
DO ii=1,ndims
    istat = NF90_INQUIRE_DIMENSION(ncid,ii,dimnames(ii),lendims(ii))
    ! Select case seems not work as we compare the select variable to another variable (invarnames)
    IF (dimnames(ii)=="lon" .OR. dimnames(ii)=="long" .OR.&
      dimnames(ii)=="longitude" .OR. dimnames(ii)=="g0_lon_3") THEN
       nLong=lendims(ii)
       inlon=dimnames(ii)
       istat = NF90_INQ_DIMID(ncid, inlon, LonDimID)
    ELSE IF (dimnames(ii)=="lat" .OR. dimnames(ii)=="g0_lat_2" .OR.&
      dimnames(ii)=="latitude" ) THEN
       nLats=lendims(ii)
       inlat=dimnames(ii)
       istat = NF90_INQ_DIMID(ncid, inlat, LatDimID)
    ELSE IF (dimnames(ii)=="time" .OR. dimnames(ii)=="initial_time0_hours") THEN
       nTime=lendims(ii)
       intime=dimnames(ii)
       istat = NF90_INQ_DIMID(ncid, intime, TimeDimID)
    ELSE IF (dimnames(ii)=="lev" .OR. dimnames(ii)== "level" .OR. dimnames(ii)=="lv_ISBL1") THEN
       nVert=lendims(ii)
       inver=dimnames(ii)
       istat = NF90_INQ_DIMID(ncid, inver, VerDimID)
       IF (istat.EQ.nf90_NoErr) exlev=.TRUE.
    ELSE IF (dimnames(ii)=="time_bnds" .OR. dimnames(ii)=="bnds") THEN
       PRINT*, "Omit time_bnds"
    ELSE 
       WRITE(*,*) '-----------------------------------------'
       WRITE(*,*) '----- Unknown dimension, check file -----'
       STOP
   ENDIF  ! IF dimnames(dimn)   
ENDDO

IF(.NOT.exlev) nVert=-1

ALLOCATE(times(nTime), latitude(nLats), longitude(nLong),vertical(nVert))
ALLOCATE(arr(nLong,nLats,tlen),arr2(nLong,nLats,tlen))
ALLOCATE(arr3(nLong,nLats,tlen))
IF (nVert.NE.-1) ALLOCATE(arrin(nLong,nLats,nVert,tlen))

istat = NF90_INQ_VARID(ncid, intime, TimeVarID)
IF(istat /= nf90_NoErr) WRITE(*,*) "Input ERROR: time dimension: inq_varid: ",NF90_STRERROR(istat)
istat = NF90_GET_VAR(ncid, TimeVarID, times)
IF(istat /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(istat)
istat = NF90_GET_ATT(ncid, TimeVarID, 'units', cftimeunit)
IF(istat /= nf90_NoErr) WRITE(*,*) 'Input ERROR: time dimension: units attribute:', NF90_STRERROR(istat)
istat = NF90_GET_ATT(ncid, TimeVarID, 'calendar', cfcalendar)
IF(istat /= nf90_NoErr) WRITE(*,*) 'Input ERROR:time dimension: calendar attribute:', NF90_STRERROR(istat)

istat = NF90_INQ_VARID(ncid, inlon, LonVarID)
IF(istat /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(istat)
istat = NF90_GET_VAR(ncid, LonVarID, longitude)
IF(istat /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(istat)

istat = NF90_INQ_VARID(ncid, inlat, LatVarID)
IF(istat /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(istat)
istat = NF90_GET_VAR(ncid, LatVarID, latitude)
IF(istat /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(istat)
IF (exlev) THEN
   istat = NF90_INQ_VARID(ncid, inver, VerVarID)
   IF(istat /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(istat)
   istat = NF90_GET_ATT(ncid, VerVarID, 'units', cflevunit)
   IF(istat /= nf90_NoErr) WRITE(*,*) 'Input ERROR:vertical dimension: units attribute:', NF90_STRERROR(istat)
   istat = NF90_GET_ATT(ncid, VerVarID, 'long_name', cflevname)
   IF(istat /= nf90_NoErr) WRITE(*,*) 'Input ERROR:vertical dimension: long_name attribute:', NF90_STRERROR(istat)
   istat = NF90_GET_VAR(ncid, VerVarID, vertical)
   IF(istat /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(istat)
ELSE
   cflevunit='NULL'
   cflevname='NULL'
   vertical=1
ENDIF

ALLOCATE(arr2d(nLong,nLats))
istat = NF90_INQ_VARID(ncid, invarname, DatVarID)
IF(istat /= nf90_NoErr) THEN 
   WRITE(*,*) NF90_STRERROR(istat)
   WRITE(*,*) "-->Most probably to second argument (invar) is wrong --> see ./tm2d -h"
   STOP
ENDIF

!==================================================================================================================================
!=========================Prepare output file


cfoutcomment='TM2D version '//trim(version)
istat=0
IF (compression) THEN
   istat = NF90_CREATE(trim(outfile), NF90_HDF5, ncido)
ELSE
   istat = NF90_CREATE(trim(outfile), NF90_NOCLOBBER, ncido)
ENDIF ! if compression [ change whether zipped NC4 file is created or generic NC3 (much larger!)
IF(istat /= nf90_NoErr) THEN ; WRITE(*,*) 'Create file: ', NF90_STRERROR(istat); ENDIF

! Define Dimension
istat = NF90_DEF_DIM(ncido,"lon",nLong,LonDimIDo)  ; IF(istat/=NF90_NoErr) print*,'Create Lon', NF90_STRERROR(istat)
istat = NF90_DEF_DIM(ncido,"lat",nLats,LatDimIDo)  ; IF(istat/=NF90_NoErr) print*,'Create Lat', NF90_STRERROR(istat)
istat = NF90_DEF_DIM(ncido,"time",nf90_unlimited,TimeDimIDo); IF(istat/=NF90_NoErr) print*,'Create Time',NF90_STRERROR(istat)

! Define Variables
istat = NF90_DEF_VAR(ncido,"lon",NF90_FLOAT,(/ LonDimIDo /),LonVarIDo)
IF(istat/=NF90_NoErr) print*,'Define lon', NF90_STRERROR(istat)
istat = NF90_PUT_ATT(ncido,LonVarIDo, "standard_name",'Longitude')
istat = NF90_PUT_ATT(ncido,LonVarIDo, "long_name",'Longitude')
istat = NF90_PUT_ATT(ncido,LonVarIDo, "units",'degree_east')


istat = NF90_DEF_VAR(ncido,"lat",NF90_FLOAT,(/ LatDimIDo /),LatVarIDo)
IF(istat/=NF90_NoErr) print*,'Define lat', NF90_STRERROR(istat)
istat = NF90_PUT_ATT(ncido,LatVarIDo, "standard_name",'Latitude')
istat = NF90_PUT_ATT(ncido,LatVarIDo, "long_name",'Latitude')
istat = NF90_PUT_ATT(ncido,LatVarIDo, "units",'degree_north')

istat = NF90_DEF_VAR(ncido,"time",NF90_DOUBLE,(/ TimeDimIDo /),timeVarIDo)
IF(istat/=NF90_NoErr) print*,'Define time ', NF90_STRERROR(istat)
istat = NF90_PUT_ATT(ncido,TimeVarIDo, "standard_name",'Time')
istat = NF90_PUT_ATT(ncido,TimeVarIDo, "long_name",'Time')
istat = NF90_PUT_ATT(ncido,TimeVarIDo, "units",cftimeunit)
istat = NF90_PUT_ATT(ncido,TimeVarIDo, "calendar",cfcalendar)
IF(istat/=NF90_NoErr) print*,'Define time ', NF90_STRERROR(istat)

! Define 3D output file
IF (labeling) THEN
   istat = NF90_DEF_VAR(ncido,"blocks",NF90_INT,(/LonDimIDo,LatDimIDo,TimeDimIDo/), DatVarIDo)
ELSE 
   istat = NF90_DEF_VAR(ncido,"blocks",NF90_BYTE,(/LonDimIDo,LatDimIDo,TimeDimIDo/), DatVarIDo)
ENDIF ! If numbering is requested, we need more precision to store the whole range of numbers (not only binary 0/1 field)

IF(istat/=NF90_NoErr) print*,'Define outdat ', NF90_STRERROR(istat)
IF (compression) THEN
   istat = NF90_DEF_VAR_DEFLATE(ncido,DatVarIDo,0,1,1)
   IF(istat/=NF90_NoErr) print*,'Define deflation ', NF90_STRERROR(istat)
ENDIF ! Compress NC4 file if compression is turned on
istat = NF90_PUT_ATT(ncido,DatVarID, "standard_name",'blocks')
istat = NF90_PUT_ATT(ncido,DatVarID, "long_name",'blocks')
istat = NF90_PUT_ATT(ncido,DatVarID, "units",'-')
istat = NF90_PUT_ATT(ncido,DatVarID, "missing_value",-127)

istat = NF90_PUT_ATT(ncido,NF90_GLOBAL, 'title', 'TM2D blocking detection algorithm')
istat = NF90_PUT_ATT(ncido,NF90_GLOBAL, 'institution','University of Bern: Group Broennimann')

istat = NF90_ENDDEF(ncido)
IF(istat/=NF90_NoErr) print*,'Error while closing definition mode: ', NF90_STRERROR(istat)

! Check whether latitude is reversed or not:
IF (latitude(1)-latitude(nLats).LT.0) THEN
   WRITE(*,'(A19,F7.3,A4,F7.3,A17)') "Latitude goes from ",latitude(1),"to ",latitude(nLats),". --> Reverse it!"
   latitude=latitude(nLats:1:-1)
   reverselat=.TRUE.
ENDIF

! Put data in NetCDF file
istat = NF90_PUT_VAR(ncido,LonVarIDo,longitude)
IF(istat/=NF90_NoErr) print*,'Put lon ', NF90_STRERROR(istat)
istat = NF90_PUT_VAR(ncido,LatVarIDo,latitude)
IF(istat/=NF90_NoErr) print*,'Put lat ', NF90_STRERROR(istat)
!istat = NF90_PUT_VAR(ncido,TimeVarID,REAL(times))
istat = NF90_PUT_VAR(ncido,TimeVarIDo,times,start=(/1/),count=(/nTime/))
IF(istat/=NF90_NoErr) print*,'Put time: ', NF90_STRERROR(istat)


!==================================================================================================================================
!==========================Allocate array for contour tracking and output
ALLOCATE(iiv(nLong*nLats*100),jjv(nLong*nLats*100))
ALLOCATE(ttv(nLong*nLats*100),dzv(nLong*nLats*100))
ALLOCATE(ttva(nLong*nLats*100))
ALLOCATE(arrover(nLong,nLats,3))
IF (logging) ALLOCATE(dzmaxv(nTime),dzmaxxv(nTime),dzmaxyv(nTime))

IF (logging) ALLOCATE(btrack(nTime*10,13)) ! Tests show that we have ~4 blocks/timestep
IF (logging) ALLOCATE(bstat(INT(nTime*1.5),12))

!==================================================================================================================================
! =========================Start calculation

!! Get dx and dy and setup tracking routine
dx=ABS(longitude(1)-longitude(2)) !360./REAL(nLong)
dy=ABS(latitude(1)-latitude(2)) !180./REAL(nLats-1)
IF (mode=='TM2D') THEN
   sy=latitude(1)
   IF(sy.GT.30) nh=.TRUE.      ! IF only SH data don't calculate NH
   ey=latitude(nLats)
   IF(ey.LT.-30) sh=.TRUE.     ! IF only NH data don't calculate SH
   nstep=INT((90-minlat)/dy)+1 ! Define number of iterations depending on resolution and most equatorward lat
   tmo=INT(tmy/dy)             ! Offset needed between the upper and the lower grid point chosen for comparison
ENDIF ! Define variables relevant for TM2D

! In the next step we load data from the NetCDF file. Because the whole file may be very large, the program may fail due to RAM 
! shortage. Therefore we only load a chunk of data (defined by tlen) into the RAM and continue to overwrite it as we move forward
! in the file. tlen should be long enough, so that no block will be longer than tlen. Such a block is not completely in the loaded
! data (i.e. in the arr variable) and the algorithm is failing.

DO tt=1,(nTime+tlen)
   IF (tt.LE.nTime) THEN 
      IF (MOD(times(tt),10000.D0)-101.00.LT.0.001) WRITE(*,'(A14,I5)') "Calculate year",INT(times(tt)/10000)
   ENDIF ! tt
   tx=MOD(tt,tlen)
   IF(tx.EQ.0) tx=tlen
   ! Load one time step slice from netcdf file, as long as we haven't reached the end of the NC file (nTime)
   IF (tt.LE.nTime) THEN
      IF (exlev) THEN
         istat = NF90_GET_VAR(ncid, DatVarID, arr2d,start = (/1,1,1,tt/),count=(/nLong,nLats,1,1/))
      ELSE
         istat = NF90_GET_VAR(ncid, DatVarID, arr2d,start = (/1,1,tt/),count=(/nLong,nLats,1/))
      ENDIF ! IF exlev   
      IF(istat /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(istat)

      IF (reverselat) arr2d=arr2d(:,nLats:1:-1) ! Reverse latitude if necessary
      arr(:,:,tx)=arr2d
   ELSE ! if end of file is reached, we load dummy data
      arr(:,:,tx)=-99999
   ENDIF ! if tt
   arr2(:,:,tx)=0; arr3(:,:,tx)=0  ! Reset arr2/arr3, otherwise we have the blockings from the former iterations in the matrix

!! ======================================================================================
!! Get potentially blocked grid points --> Depends on the applied method --> mode switch
   IF (tt.LE.nTime) THEN ! We do not need to calculate blockings after the last time step
      ! Get potentially blocked grid points for TM2D
      IF (mode=="TM2D") THEN
         IF (nh) THEN
            IF(tt.EQ.1) print*,"We do NH --> Tibaldi 1990 (changed to 15°)"
            DO ii=1,nLong
               DO jj=tmo+1,nstep
                  !    Zn           -    Zo         < -10*deg (see TM90)
                  IF(arr(ii,(jj-tmo),tx)-arr(ii,jj,tx).LT.(tmy*dzp).AND.& ! GHGN (polar)
                  !    Zo           -    Zs         > 0*deg
                  arr(ii,jj,tx)-arr(ii,(jj+tmo),tx).GT.(tmy*dze)) THEN    ! GHGS (equator)
                     arr2(ii,jj,tx)=1
                  ENDIF
               ENDDO ! jj
            ENDDO ! ii
         ENDIF ! IF nh

         IF (sh) THEN
            IF(tt.EQ.1) print*,"We do SH --> Tibaldi 1994"
            DO ii=1,nLong !nLats-nstep-tmo,nLats-tmo
               DO jj=(nLats-nstep+1),(nLats-tmo)
                  !         Zn(polar)    - Zo        < -10*deg           
                  IF(arr(ii,(jj+tmo),tx)-arr(ii,jj,tx).LT.(tmy*dzp).AND.&
                  !         Zo    - Zs (equator)     > 0*deg
                  arr(ii,jj,tx)-arr(ii,(jj-tmo),tx).GT.(tmy*dze)) THEN
                     arr2(ii,jj,tx)=1
                  ENDIF
               ENDDO ! jj
            ENDDO ! ii
         ENDIF ! IF sh
      ! Get potentially blocked grid points for Z500anom
      ELSE IF (mode=="Z500anom") THEN
         IF(tt.EQ.1) print*,"Get anomalies"
         DO ii=1,nLong
            DO jj=1,nLats
               IF(arr(ii,jj,tx).GT.anom) THEN
                  arr2(ii,jj,tx)=1
               ENDIF ! if PV < -1.3
            ENDDO ! jj
         ENDDO ! ii
      ! Get potentially blocked grid points for VAPV
      ELSE IF (mode=="VAPV") THEN
         IF(tt.EQ.1) print*,"Get PV anomalies"
         DO ii=1,nLong
            DO jj=1,nLats
               IF(latitude(jj).LT.0) THEN
                  arr(ii,jj,tx)=arr(ii,jj,tx)*(-1)
               ENDIF ! IF southern hemisphee, change sign of PV
               IF(arr(ii,jj,tx).LT.vapv) THEN
                   arr2(ii,jj,tx)=1
               ENDIF ! if PV < -1.3
            ENDDO ! jj
         ENDDO ! ii
      ELSE
         PRINT*,"ERROR: mode not found --> Aborting"
         STOP
      ENDIF ! if mode

!! Identify individual contour areas on each lon-lat field
      ALLOCATE(carr(nLong,nLats))
      IF(tt.EQ.1) print*,"Contour2D"
      carr=arr2(:,:,tx)
      nclu=0
      label2d=1
      DO ii=1,nLong
         DO jj=1,nLats
            ! Find a potentially blocked grid point and connect all spatially connected blocked grid points
            IF (carr(ii,jj).EQ.1) THEN
               label2d=label2d+1
               carr(ii,jj)=label2d
               CALL connect2D(carr,label2d,ii,jj,nLong,nLats)
            ENDIF ! IF carr
         ENDDO ! jj
      ENDDO ! ii
      arr2(:,:,tx)=carr
      DEALLOCATE(carr)
   ENDIF ! IF arr(1,1,tx)

!! Calculate overlap and remove items where overlap is too small
   IF(tt.EQ.1) print*,"Overlap"
   t1=tx-2 ; t2=tx-1 ; t3=tx
   IF(t1.LT.1) t1=t1+tlen
   IF(t2.LT.1) t2=t2+tlen
   IF(t3.GT.tlen) t3=t3-tlen

   IF(tt.EQ.2) THEN 
      arrover(:,:,1)=0 ; arrover(:,:,2:3)=arr2(:,:,(/t2,t3/))
   ELSEIF (tt.EQ.(nTime+1)) THEN
      arrover(:,:,3)=0 ; arrover(:,:,1:2)=arr2(:,:,(/t1,t2/)) 
   ELSEIF (tt.GT.2.AND.tt.LE.nTime) THEN
      arrover=arr2(:,:,(/t1,t2,t3/))
   ENDIF ! Get array for overlap

   IF (tt.GT.1.AND.tt.LT.(nTime+1)) THEN 
      ! Here we calculate the overlap --> overlap defined as overlapping area of t2 and t3 compared to total area in t2
      CALL spoverlap(arrover,nLong,nLats,latitude,dx,dy,overlap)
      arr2(:,:,t2)=arrover(:,:,2)
      arr3(:,:,t2)=arrover(:,:,2)
   ENDIF ! Define overlap for every timestep --> Load 3 timesteps in spoverlap and compute ovelap of the middle timestep to the first and third


!! Identify individual contour areas in 3D field (lon-lat-time)
   IF(tt.EQ.1) print*,"Contour3D"
   IF (tt.GT.1.AND.tt.LE.(nTime+1)) THEN 
      WHERE(arr2(:,:,t2).GE.1) arr3(:,:,t2)=-1
   ENDIF ! IF tt.GT.1

   IF(tt.GT.(tlen-10).AND.tt.LE.nTime+(tlen-10)) THEN
      t1=tt-(tlen-10)
      t2=tt-(tlen-10) 
      IF (t1.LT.1) t1=t1+tlen
      DO WHILE (t2.GT.tlen) 
         t2=t2-tlen
      ENDDO ! DO WHILE
      DO ii=1,nLong
         DO jj=1,nLats
            IF (arr3(ii,jj,t2).EQ.-1) THEN
               label=label+1
               arr3(ii,jj,t2)=label
               ngpc=1       ; ngp=1
               IF (logging) THEN 
                 dzmax=0      ; dzmaxv=0  ; dzmaxxv=0 ; dzmaxyv=0
               ENDIF ! LOGGING
               ! Get values for first detected grid point
               iiv(1)=ii    ; jjv(1)=jj ; ttv(1)=t2 ; ttva(1)=t1
               IF (logging) THEN
                  IF (mode=="TM2D") THEN 
                     IF (jj.LT.(nLats/2)) dzv(1)=arr(ii,jj,t2)-arr(ii,(jj+tmo),t2) ! NH
                     IF (jj.GE.(nLats/2)) dzv(1)=arr(ii,jj,t2)-arr(ii,(jj-tmo),t2) ! SH
                  ELSE IF ((mode=="VAPV") .OR. (mode=="Z500anom")) THEN
                     dzv(1)=arr(ii,jj,t2)
                  ENDIF ! IF mode
                  ! Get data for blockstat
                  dzmax=dzv(1) ; dzmaxt=t1 ; dzmaxx=ii ; dzmaxy=jj
                  ! Get data for blocktracks
                  dzmaxv(t1)=dzv(1) ; dzmaxxv(t1)=ii ; dzmaxyv(t1)=jj
               ENDIF ! IF logging

               DO WHILE (ngpc.LE.ngp) ! Search blocked points around each point. Add new points to stack until we're done with all blocked points
                  DO ik=-1,1
                     DO jl=-1,1
                        DO tm=-1,1
                           ! Conditions to inhibit to go out of bounds of array
                           IF(ABS(ik)+ABS(jl)+ABS(tm).NE.1) CYCLE
                           ix=iiv(ngpc)+ik
                           IF(ix.EQ.0)       ix=nLong
                           IF(ix.GT.nLong)   ix=1
                           iy=jjv(ngpc)+jl
                           IF(iy.EQ.0.OR.iy.GT.nLats) CYCLE
                           it=ttv(ngpc)+tm
                           iu=ttva(ngpc)+tm
                           IF (iu.LT.ttv(1)) CYCLE ! We don't want to go backward where step t+(nlen-1) is already loaded (absolute t-axis!)
                           IF (it.GT.tlen)   it=it-tlen ! If we reach the end of the array, go to first entry
                           IF (it.LT.1)      it=it+tlen
                 
                           IF (arr3(ix,iy,it).EQ.-1) THEN ! Add point to stack if not already in stack
                              arr3(ix,iy,it)=label
                              ngp=ngp+1
                              iiv(ngp)=ix ; jjv(ngp)=iy ; ttv(ngp)=it ; ttva(ngp)=iu
                              
                              IF (logging) THEN
                                 ! Get maximum anomaly depending on the chose mode
                                 IF (mode=="TM2D") THEN
                                    IF (jj.LT.(nLats/2)) dzv(ngp)=arr(ix,iy,it)-arr(ix,(iy+tmo),it) ! NH
                                    IF (jj.GE.(nLats/2)) dzv(ngp)=arr(ix,iy,it)-arr(ix,(iy-tmo),it) ! SH
                                 ELSE IF ((mode=="Z500anom") .OR. (mode=="VAPV")) THEN
                                    dzv(ngp)=arr(ix,iy,it)
                                 ENDIF

                                 ! Get data for blockstata
                                 IF ((mode=="TM2D") .OR. (mode=="Z500anom")) THEN
                                    ! Get data for blockstat
                                    IF (dzv(ngp).GT.dzmax) THEN
                                       dzmax=dzv(ngp)
                                       dzmaxt=iu
                                       dzmaxx=ix
                                       dzmaxy=iy
                                    ENDIF ! Dzmin
                                    ! Get data for blocktracks
                                    IF(dzv(ngp).GT.dzmaxv(iu)) THEN
                                       dzmaxv(iu)=dzv(ngp)
                                       dzmaxxv(iu)=ix
                                       dzmaxyv(iu)=iy
                                    ENDIF ! dzmaxv
                                 ELSE IF (mode=="VAPV") THEN
                                    ! Get data for blockstat
                                    IF (dzv(ngp).LT.dzmax) THEN ! Changed for PV ( .GT: for GPH
                                       dzmax=dzv(ngp)
                                       dzmaxt=iu
                                       dzmaxx=ix
                                       dzmaxy=iy
                                    ENDIF ! Dzmin
                                    ! Get data for blocktracks
                                    IF(dzv(ngp).LT.dzmaxv(iu)) THEN ! Changed for PV ( .GT. for GPH)
                                       dzmaxv(iu)=dzv(ngp)
                                       dzmaxxv(iu)=ix
                                       dzmaxyv(iu)=iy
                                    ENDIF ! dzmaxv
                                 ENDIF ! mode
                              ENDIF ! IF logging
                           ENDIF ! arr2=-1
                        ENDDO ! tm
                     ENDDO ! jl
                  ENDDO ! ik
                  ngpc = ngpc+1
               ENDDO ! WHILE

               ! Get output for blockstat
               IF (logging) THEN
                  bstat(label,1)=label                                   ! Label (Block ID) 
                  bstat(label,3)=MINVAL(ttva(1:ngp))                     ! Get start time step of block
                  bstat(label,4)=MAXVAL(ttva(1:ngp))                     ! Get end time step of block
                  bstat(label,2)=bstat(label,4)-bstat(label,3)+1         ! Get length of block
                  IF ((mode=="TM2D") .OR. (mode=="Z500anom")) THEN
                     bstat(label,5)=MAXVAL(dzmaxv)                       ! Maximum anomaly for TM2D, Z500anom
                  ELSE IF (mode=="VAPV") THEN
                     bstat(label,5)=MINVAL(dzmaxv)                       ! Minimum anomaly for PV (the lower PV, the stronger the block)
                  ENDIF ! Get maximum anomaly for block (depends on method)
                  bstat(label,6)=dzmaxt                                  ! Get time step of maximum anomaly of block
                  bstat(label,7)=dzmaxx                                  ! Get longitude index for maximum anomaly of block
                  bstat(label,8)=dzmaxy                                  ! Get latitude index for maximum anomaly of block
                  CALL gettime(times,MINVAL(ttva(1:ngp)),cftimeunit,year,mon,day,hour)
                  bstat(label,9:12)=(/year,mon,day,hour/)                ! Get year,mon,day,hour of first time step of block

                  ! Get output for blocktracks
                  ALLOCATE(sinsum(INT(bstat(label,2))),cossum(INT(bstat(label,2))))  ! bstat(label,2) contains the length of the block
                  ALLOCATE(latsum(INT(bstat(label,2))),ngpsum(INT(bstat(label,2))))
                  sinsum=0 ; cossum=0 ; latsum=0 ; ngpsum=0               ! sinsum and cossum are needed for the weighted mean of the longitude
                  DO jl=1,ngp
                     tm=ttva(jl)-MINVAL(ttva(1:ngp))+1                    ! We use the absolute time axis for calculations concerning time
                     btrack(nn+tm,11)=btrack(nn+tm,11)+REAL(ABS(12321.*dx*dy*cos(pi/180*latitude(jjv(ngp))))) ! Get area of blocking (sum over gridpoints)
                     sinsum(tm)=sinsum(tm)+sin(pi/180*longitude(iiv(jl))) ! Get sin from each gridpoint longitude
                     cossum(tm)=cossum(tm)+cos(pi/180*longitude(iiv(jl))) ! Get cos from each gridpoint longitude
                     latsum(tm)=latsum(tm)+latitude(jjv(jl))              ! Get latitude
                     ngpsum(tm)=ngpsum(tm)+1
                  ENDDO

                  DO tm=1,INT(bstat(label,2))
                     btrack(nn+tm,1)=label                                            ! Block ID / label
                     btrack(nn+tm,2)=bstat(label,2)                                   ! Blocking length
                     btrack(nn+tm,3)=tm                                               ! Time step of block
                     btrack(nn+tm,4)=dzmaxv(MINVAL(ttva(1:ngp))+tm-1)                 ! Get maximum anomaly for time step (MINVAL(ttva.. provides the time step of the first time step of block
                     CALL gettime(times,MINVAL(ttva(1:ngp)+tm-1),cftimeunit,year,mon,day,hour)
                     btrack(nn+tm,5:8)=(/year,mon,day,hour/)                          ! Get year,mon,day,hour of time step
                     btrack(nn+tm,9)=dzmaxxv(MINVAL(ttva(1:ngp))+tm-1)                ! Get longitude of maximum anomaly for time step
                     btrack(nn+tm,10)=dzmaxyv(MINVAL(ttva(1:ngp))+tm-1)               ! Get latitude of maximum anomaly for time step
                     btrack(nn+tm,12)=REAL(ATAN2(sinsum(tm),cossum(tm))/pi*180)       ! Calculate (longitude) angle out of the sine and cosine component
                     IF (btrack(nn+tm,12).LT.0) btrack(nn+tm,12)=btrack(nn+tm,12)+360 ! If negative, make a positive number
                     btrack(nn+tm,13)=latsum(tm)/ngpsum(tm) ! Get mean latitude       ! Calculate mean (latitude) for time step
                  ENDDO
                  DEALLOCATE(sinsum,cossum,latsum,ngpsum)

                  nn=nn+INT(bstat(label,2))
               ENDIF ! logging

               ! Persistence criterium - delete labels out of array if block is too short
               IF ((MAXVAL(ttva(1:ngp))-MINVAL(ttva(1:ngp))+1).LT.persistence) THEN
                  DO tm=1,ngp
                     arr3(iiv(tm),jjv(tm),ttv(tm))=0
                     arr2(iiv(tm),jjv(tm),ttv(tm))=0
                  ENDDO ! DO erase all blocked points that are not persistent  
                  !label=label-1
               ENDIF ! If persistence
               ttv=0; iiv=0; jjv=0

            ENDIF ! arr3.EQ.-1
         ENDDO ! DO jj (lats)
      ENDDO ! DO ii (lons)

      IF ( .NOT. labeling) THEN
          WHERE(arr3(:,:,t2).GT.1) arr3(:,:,t2)=1
      ENDIF ! Set blocked areas to 1, in case blocks numbering/labeling is not desired (standard)
      arr2d=arr3(:,:,t2)

      ! Put time slice into NetCDF file
      istat = NF90_PUT_VAR(ncido,DatVarIDo,arr2d,start = (/1,1,t1/),count=(/nLong,nLats,1/))
   ENDIF ! IF(tt.GT.(tlen-10)) THEN


ENDDO ! DO tt
istat = NF90_CLOSE(ncid)
istat = NF90_CLOSE(ncido)
print*,"logging",logging
IF (logging) THEN 
   OPEN(21,FILE="blockstat_all.txt",ACTION="write",STATUS="unknown",POSITION="rewind")
   OPEN(22,FILE="blockstat.txt",ACTION="write",STATUS="unknown",POSITION="rewind")
   OPEN(23,FILE="blocktracks.txt",ACTION="write",STATUS="unknown",POSITION="rewind")
   OPEN(24,FILE="blocktracks_all.txt",ACTION="write",STATUS="unknown",POSITION="rewind")


   900 FORMAT (2(A8),A5,3(A3),2(A8),A10,  A8,2(A5),2(A9))
   901 FORMAT (2(I8),I5,3(I3),2(I8),F10.3,I8,2(I5),2(F9.3))
   902 FORMAT (A8,2(A7),A10,  A5,3(A3),2(A9),  A11,2(A9))
   903 FORMAT (I8,2(I7),F10.3,I5,3(I3),2(F9.3),I11,2(F9.3))

   IF (mode=="TM2D") THEN
      WRITE(21,900) 'ID',"Length","YYYY","MM","DD","HH","Start","End","Max_dz","Timx","Lonx","Latx","Lon","Lat"
      WRITE(22,900) 'ID',"Length","YYYY","MM","DD","HH","Start","End","Max_dz","Timx","Lonx","Latx","Lon","Lat"
      WRITE(23,902) 'ID',"Length","Tmstp","Max_dz","YYYY","MM","DD","HH","Lonx","Latx","Area","MeanLon","MeanLat"
      WRITE(24,902) 'ID',"Length","Tmstp","Max_dz","YYYY","MM","DD","HH","Lonx","Latx","Area","MeanLon","MeanLat"
   ELSE IF (mode=="Z500anom") THEN
      WRITE(21,900) 'ID',"Length","YYYY","MM","DD","HH","Start","End","Max_anom","Timx","Lonx","Latx","Lon","Lat"
      WRITE(22,900) 'ID',"Length","YYYY","MM","DD","HH","Start","End","Max_anom","Timx","Lonx","Latx","Lon","Lat"
      WRITE(23,902) 'ID',"Length","Tmstp","Max_anom","YYYY","MM","DD","HH","Lonx","Latx","Area","MeanLon","MeanLat"
      WRITE(24,902) 'ID',"Length","Tmstp","Max_anom","YYYY","MM","DD","HH","Lonx","Latx","Area","MeanLon","MeanLat"
   ELSE IF (mode=="VAPV") THEN
      WRITE(21,900) 'ID',"Length","YYYY","MM","DD","HH","Start","End","Min_PV","Timx","Lonx","Latx","Lon","Lat"
      WRITE(22,900) 'ID',"Length","YYYY","MM","DD","HH","Start","End","Min_PV","Timx","Lonx","Latx","Lon","Lat"
      WRITE(23,902) 'ID',"Length","Tmstp","Min_PV","YYYY","MM","DD","HH","Lonx","Latx","Area","MeanLon","MeanLat"
      WRITE(24,902) 'ID',"Length","Tmstp","Min_PV","YYYY","MM","DD","HH","Lonx","Latx","Area","MeanLon","MeanLat"
   ENDIF ! IF mode
   nn=0
   DO ii=1,label
      WRITE(21,901) INT(bstat(ii,(/1,2,9,10,11,12,3,4/))),bstat(ii,5),INT(bstat(ii,6:8)),&
                longitude(INT(bstat(ii,7))),latitude(INT(bstat(ii,8))) 
      DO tm=1,INT(bstat(ii,2))
         WRITE(24,903) INT(btrack(nn+tm,1:3)),btrack(nn+tm,4),INT(btrack(nn+tm,5:8)),&
                   longitude(INT(btrack(nn+tm,9))),latitude(INT(btrack(nn+tm,10))),INT(btrack(nn+tm,11)),btrack(nn+tm,12:13)
      ENDDO ! tm 
      IF (bstat(ii,2).GE.persistence) THEN
         write(22,901) INT(bstat(ii,(/1,2,9,10,11,12,3,4/))),bstat(ii,5),INT(bstat(ii,6:8)),&
                   longitude(INT(bstat(ii,7))),latitude(INT(bstat(ii,8)))
         DO tm=1,INT(bstat(ii,2))
            IF(INT(btrack(nn+tm,6)).EQ.0) btrack(nn+tm,6)=1
            IF(INT(btrack(nn+tm,7)).EQ.0) btrack(nn+tm,7)=1 
            write(23,903) INT(btrack(nn+tm,1:3)),btrack(nn+tm,4),INT(btrack(nn+tm,5:8)),&
                           longitude(INT(btrack(nn+tm,9))),latitude(INT(btrack(nn+tm,10))),INT(btrack(nn+tm,11)),btrack(nn+tm,12:13)
         ENDDO ! tm 
      ENDIF ! IF length > persistencea
      nn=nn+INT(bstat(ii,2))
   ENDDO
   CLOSE(21)
   CLOSE(22)
   CLOSE(23)
   CLOSE(24)
   DEALLOCATE(arr,bstat,btrack)
ENDIF ! IF logging 
DEALLOCATE(iiv,jjv,ttv)

DEALLOCATE(arr2,times,longitude,latitude)
END PROGRAM tm2d

!======================== Subroutines =================================================
SUBROUTINE gettime(times,tt,cftimeunit,year,mon,day,hour)
   IMPLICIT NONE

   INTEGER, INTENT(IN)                           :: tt
   REAL(8), DIMENSION(:), INTENT(INOUT)          :: times
   CHARACTER(1000),INTENT(IN)                    :: cftimeunit
   INTEGER,INTENT(OUT)                           :: year,mon,day,hour
   
   IF (cftimeunit(1:16).EQ."day as %Y%m%d.%f") THEN
      year=INT(times(tt)/10000)
      mon=INT(MOD(INT(times(tt)),10000)/100)
      day=MOD(INT(times(tt)),100)
      hour=INT(MOD(times(tt),1.D0)*24)
   ELSE
      IF (tt.EQ.1) WRITE(*,*) "WARNING: Use an absolute time axis,the timestamp in the log files will be wrong&
                               &--> Consider preprocessing the NetCDF with cdo -a copy ifile ofile"
   ENDIF
END SUBROUTINE gettime


SUBROUTINE connect2D(carr,label,ii,jj,nLong,nLats)
   IMPLICIT NONE

   INTEGER, INTENT(IN)                              :: label, ii,jj,nLong,nLats
   INTEGER, INTENT(INOUT), DIMENSION(nLong,nLats)   :: carr

   INTEGER, DIMENSION(nLong*nLats)                  :: iiv, jjv                       ! Get locations of points in 
   INTEGER                                          :: ngp                            ! Number of gridpoints belonging to 
   INTEGER                                          :: ngpc                           ! Number of grid p
   INTEGER                                          :: ix,iy
   INTEGER                                          :: ik,jl


   ngpc=1                ! We get the coordinates for the first grid point
   iiv(1)=ii
   jjv(1)=jj
   ngp=1                 ! At the moment we have just one grid point in the contour

   ! We search for unassigned gridpoint around the initial gridpoint, if we add them to the stack (locations in iiv,jjv)
   ! After we found (several) neighbouring gridpoints which also belong to the contour, we check for neighbouring grid points
   ! for each of these gridpoints. As soon as we find no new gridpoints, the DO WHILE loop stops
   DO WHILE ( ngpc.LE.ngp)
      DO ik=-1,1
         DO jl=-1,1
            IF(ABS(ik)+ABS(jl).GT.1) CYCLE
            ix=iiv(ngpc)+ik 
            IF(ix.EQ.0) THEN 
               ix=nLong 
            ELSEIF (ix.GT.nLong) THEN 
               ix=ix-nLong 
            ENDIF ! IF ix --> Correct ix (longitude) when we are at the edge of the matrix
            iy=jjv(ngpc)+jl ; IF(iy.EQ.0.OR.iy.GT.nLats) CYCLE ! IF iy --> Inhibit iy (latitude) to get out of bounds at N/S pole

            IF ( (carr(ix,iy).EQ.1).AND.(carr(ix,iy).NE.label) ) THEN
               carr(ix,iy)=label
               ngp        = ngp+1 
               iiv(ngp)   = ix !iiv(ngpc)+ik ! Save lon of newly found grip point
               jjv(ngp)   = iy !jjv(ngpc)+jl ! Save lat of newly found grip point
            ENDIF ! IF gridpoint not yet included
         ENDDO ! jj
      ENDDO ! ii
      ngpc = ngpc+1 ! ngpc gets larger than ngp as soon as we don't add new gridpoints
   ENDDO ! DO while ngpc <= ngp
END SUBROUTINE connect2D

SUBROUTINE spoverlap(arrover,nLong,nLats,latitude,dx,dy,overlap)
! Compute spatial overlap between time step 2 and 3 and compared it to area of block at time step 2
! Time step 1 is currently unused but remains in the code in order to make changes to the overlapping algorithm easier

IMPLICIT NONE

   INTEGER,INTENT(IN)                                   :: nLong,nLats
   REAL,INTENT(IN)                                      :: dx,dy,overlap
   REAL,INTENT(IN),DIMENSION(nLats)                     :: latitude
   INTEGER,INTENT(INOUT),DIMENSION(nLong,nLats,3)       :: arrover
   INTEGER                                              :: nn,ii,jj
   REAL(8)                                              :: pi=4.0*ATAN(1.D0)
   REAL(8)                                              :: a1,ax,afrc

   ! Loop over every labelled block of time step 2
   DO nn=2,INT(MAXVAL(arrover(:,:,2))) ! For each contour
      a1=0 ; ax=0 ; afrc=0 
      DO ii=1,nLong
         DO jj=1,nLats
            ! Get area of contour at time tt and tt+1 --> Overlap between both timesteps (forward)

            ! Get total area of time step 2
            IF (arrover(ii,jj,2).EQ.nn) a1=a1+ABS(12321.*dx*dy*cos(pi/180*latitude(jj)))
            ! Get area that is overlapping with time step 3
            IF ((arrover(ii,jj,2).EQ.nn).AND.(arrover(ii,jj,3).GE.1)) ax=ax+ABS(12321.*dx*dy*cos(pi/180*latitude(jj)))
            afrc=ax/a1
         ENDDO ! jj
      ENDDO !! ii

      IF (afrc.LT.overlap) THEN
         WHERE(arrover(:,:,2).EQ.nn) arrover(:,:,2)=0
      ENDIF ! IF overlap not big enough
   ENDDO ! nn
END SUBROUTINE spoverlap


