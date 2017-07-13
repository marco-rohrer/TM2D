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


! Marco Rohrer, Oct 2016 (marco.rohrer@giub.unibe.ch)


USE netcdf

IMPLICIT NONE


LOGICAL                                    :: nopersistence=.FALSE. !.TRUE. !FALSE.
LOGICAL                                    :: nooverlap=.FALSE. !.TRUE. !.FALSE.

! NetCDF I/O
CHARACTER(1000)                            :: cfoutcomment, cftimeunit, cfcalendar             ! Variables to store NetCDF comments and time units
CHARACTER(200)                             :: infile, invarname, inlat, inlon, inver, intime   ! Names of the incoming file, the variable name that should be read and the dimensions
CHARACTER(200)                             :: outfile                                          ! Name of outgoing file
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


!logging=.FALSE.
version='v0.8.6'
IF(iargc()==1) THEN
   CALL getarg(1,infile)
   SELECT CASE (infile)
      CASE("-h","--help","-H","--Help")
         PRINT*,"Usage tm2d infile invar outfile"
         PRINT*,"Example: ./tm2d eraint.nc z blocks_eraint.nc"
         PRINT*," "
         WRITE(*,'(A257,I3,A50)') "This programm detects and tracks blockings. First geopotential height reversal are detected wi&
                &th the same approach as Tibaldi and Molteni (1990). See also Tibaldi et al. (1994) for Southern Hemisphere. The &
                &latitudes for the calculation of the gradients are ",tmy," apart. This can be changed via the parameter tmy."
         WRITE(*,'(A142,F3.1,A115)') "Afterwards the overlap is between two time steps is used as the stationarity criterion. If &
                &the overlap between two timesteps must be at least ",overlap,",otherwise the anomaly is not counted as a consecu&
                &tive block. The overlap can be changed via the parameter overlap."
         WRITE(*,'(A175,I3,A122)') "The blockings are then detected as a spatially and temporally connected blob of grid points. &
                &If the temporal extent is larger than the minimum persistence criterium (at least, ",persistence," time steps),&
                & the blob is detected as a blocking. The persistence criterium can be changed with the persistence criterium."
         WRITE(*,'(A345)') "The output contains the NetCDF file with the name given as the third argument and four files containi&
                &ng the blocking statistics (for all blocks and for blocks reaching the persistence criterium) and the block trac&
                &ks which indicates the position and strenght of each blocking at each time step (again for all blocks and for bl&
                &ocks > persistence)."
         WRITE(*,'(A70)') "Writing logfiles can be disabled by setting the seventh parameter to F"
         WRITE(*,'(A139)') "Block tracks has two measures of the blocking location: The location of the maximum anomaly and the l&
                &ocation of the center of the blocking."
         PRINT*,""
         WRITE(*,'(A135)') "Optionally values all three tuning variables can be given after the first three arguments in the foll&
                &owing form: persistence overlap latdiff"  
         WRITE(*,'(A58)') "==> Example: ./tm2d eraint.nc z blocks_eraint.nc 20 0.7 15"
         WRITE(*,'(A70)') "Additionally the logfiles can be turned off by giving a 7th argument F:"
         WRITE(*,'(A60)') "==> Example: ./tm2d eraint.nc z blocks_eraint.nc 20 0.7 15 F"
      CASE("--Version","-V","--version")
            PRINT*, "Version: ",version
            PRINT*, "Written by M. Rohrer (marco.rohrer@giub.unibe.ch),last changed 20161007"
    END SELECT
    STOP
ENDIF

IF(iargc().LT.3.OR.iargc().GT.7.OR.iargc().EQ.4.OR.iargc().EQ.5) THEN
    PRINT*,'Error: Wrong number of Input Arguments given. Correct form is:'
    PRINT*,'ifile.nc invar ofile.nc'
    PRINT*,'or'
    PRINT*,'ifile.nc invar ofile.nc persistence overlap latdiff'
    PRINT*,'Example: ./tm2d in.nc invar out.nc'
    PRINT*,'See ./tm2d -h for more information'
    STOP
ENDIF ! IF Iargc
  
! READ ARGUMENTS
CALL getarg(1,infile)
CALL getarg(2,invarname)
CALL getarg(3,outfile)
IF (iargc().GE.6) THEN
   CALL getarg(4,cpers)    ; READ (cpers,'(I5)')   persistence
   CALL getarg(5,cover)    ; READ (cover,'(F5.3)') overlap
   CALL getarg(6,ctmy)     ; READ (ctmy,'(I5)')  tmy
   IF (iargc().EQ.7) THEN
      CALL getarg(7,clogf) ; IF (clogf.EQ."F") logging=.FALSE.
   ENDIF ! if 7 arguments
ENDIF

! Check whether overlap criterium shall be used
IF (nooverlap) THEN
   overlap=-1
   WRITE(*,*) "ATTENTION: Overlap criterium is disabled!"
ENDIF ! nooverlap

IF (nopersistence) THEN
   persistence=0
   WRITE(*,*) "ATTENTION: Persistence criterium is disabled!"
ENDIF ! nopersistence  

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
istat = NF90_CREATE(trim(outfile), NF90_CLOBBER, ncido)
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
istat = NF90_DEF_VAR(ncido,"blocks",NF90_BYTE,(/LonDimIDo,LatDimIDo,TimeDimIDo/), DatVarIDo)
IF(istat/=NF90_NoErr) print*,'Define outdat ', NF90_STRERROR(istat)
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
sy=latitude(1)
IF(sy.GT.30) nh=.TRUE.      ! IF only SH data don't calculate NH
ey=latitude(nLats)
IF(ey.LT.-30) sh=.TRUE.     ! IF only NH data don't calculate SH
nstep=INT((90-minlat)/dy)+1 ! Define number of iterations depending on resolution and most equatorward lat
tmo=INT(tmy/dy)             ! Offset needed between the upper and the lower grid point chosen for comparison

DO tt=1,(nTime+tlen)
   IF (tt.LE.nTime) THEN 
      IF (MOD(times(tt),10000.D0)-101.00.LT.0.001) WRITE(*,'(A14,I5)') "Calculate year",INT(times(tt)/10000)
   ENDIF ! tt
   tx=MOD(tt,tlen)
   IF(tx.EQ.0) tx=tlen
   IF (tt.LE.nTime) THEN
      IF (exlev) THEN
         istat = NF90_GET_VAR(ncid, DatVarID, arr2d,start = (/1,1,1,tt/),count=(/nLong,nLats,1,1/))
      ELSE
         istat = NF90_GET_VAR(ncid, DatVarID, arr2d,start = (/1,1,tt/),count=(/nLong,nLats,1/))
      ENDIF ! IF exlev   
      IF(istat /= nf90_NoErr) WRITE(*,*) NF90_STRERROR(istat)

      IF (reverselat) arr2d=arr2d(:,nLats:1:-1) ! Reverse latitude if necessary
      arr(:,:,tx)=arr2d
   ELSE
      arr(:,:,tx)=-99999
   ENDIF ! if tt
   arr2(:,:,tx)=0; arr3(:,:,tx)=0  ! Reset arr2/arr3, otherwise we have the blockings from the former iterations in the matrix
   !! Get initial gradient inversals
!! ======================================================================================
   IF (tt.LE.nTime) THEN !(INT(arr(1,1,tx)).NE.-99999) THEN ! We do not need to calculate blockings after the last time step
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


!! Identify individual contour areas on each lon-lat field
      ALLOCATE(carr(nLong,nLats))
      IF(tt.EQ.1) print*,"Contour2D"
      carr=arr2(:,:,tx)
      nclu=0
      label2d=1
      DO ii=1,nLong
         DO jj=1,nLats
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
      CALL spoverlapnew(arrover,nLong,nLats,latitude,dx,dy,overlap)
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
                  IF (jj.LT.(nLats/2)) dzv(1)=arr(ii,jj,t2)-arr(ii,(jj+tmo),t2) ! NH
                  IF (jj.GE.(nLats/2)) dzv(1)=arr(ii,jj,t2)-arr(ii,(jj-tmo),t2) ! SH
               
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
                                 IF (jj.LT.(nLats/2)) dzv(ngp)=arr(ix,iy,it)-arr(ix,(iy+tmo),it) ! NH
                                 IF (jj.GE.(nLats/2)) dzv(ngp)=arr(ix,iy,it)-arr(ix,(iy-tmo),it) ! SH

                                 ! Get data for blockstata
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
                              ENDIF ! IF logging
                           ENDIF ! arr2=-1
                        ENDDO ! tm
                     ENDDO ! jl
                  ENDDO ! ik
                  ngpc = ngpc+1
               ENDDO ! WHILE

               ! Get output for blockstat
               IF (logging) THEN
                  bstat(label,1)=label
                  bstat(label,3)=MINVAL(ttva(1:ngp))
                  bstat(label,4)=MAXVAL(ttva(1:ngp))
                  bstat(label,2)=bstat(label,4)-bstat(label,3)+1
                  bstat(label,5)=MAXVAL(dzmaxv)
                  bstat(label,6)=dzmaxt
                  bstat(label,7)=dzmaxx
                  bstat(label,8)=dzmaxy
                  CALL gettime(times,MINVAL(ttva(1:ngp)),cftimeunit,year,mon,day,hour)
                  bstat(label,9:12)=(/year,mon,day,hour/)

                  ! Get output for blocktracks
                  ALLOCATE(sinsum(INT(bstat(label,2))),cossum(INT(bstat(label,2))))
                  ALLOCATE(latsum(INT(bstat(label,2))),ngpsum(INT(bstat(label,2))))
                  sinsum=0 ; cossum=0 ; latsum=0 ; ngpsum=0
                  DO jl=1,ngp
                     tm=ttva(jl)-MINVAL(ttva(1:ngp))+1                    ! We use the absolute time axis for calculations concerning time
                     btrack(nn+tm,11)=btrack(nn+tm,11)+REAL(ABS(12321.*dx*dy*cos(pi/180*latitude(jjv(ngp))))) ! Get area of blocking (sum over gridpoints)
                     sinsum(tm)=sinsum(tm)+sin(pi/180*longitude(iiv(jl))) ! Get sin from each gridpoint longitude
                     cossum(tm)=cossum(tm)+cos(pi/180*longitude(iiv(jl))) ! Get cos from each gridpoint longitude
                     latsum(tm)=latsum(tm)+latitude(jjv(jl))              ! Get latitude
                     ngpsum(tm)=ngpsum(tm)+1
                  ENDDO

                  DO tm=1,INT(bstat(label,2))
                     btrack(nn+tm,1)=label
                     btrack(nn+tm,2)=bstat(label,2)
                     btrack(nn+tm,3)=tm
                     btrack(nn+tm,4)=dzmaxv(MINVAL(ttva(1:ngp))+tm-1)
                     CALL gettime(times,MINVAL(ttva(1:ngp)+tm-1),cftimeunit,year,mon,day,hour)
                     btrack(nn+tm,5:8)=(/year,mon,day,hour/) 
                     btrack(nn+tm,9)=dzmaxxv(MINVAL(ttva(1:ngp))+tm-1)
                     btrack(nn+tm,10)=dzmaxyv(MINVAL(ttva(1:ngp))+tm-1)
                     btrack(nn+tm,12)=REAL(ATAN2(sinsum(tm),cossum(tm))/pi*180)  ! Calculate angle out of the sine and cosine component
                     IF (btrack(nn+tm,12).LT.0) btrack(nn+tm,12)=btrack(nn+tm,12)+360 ! If negative, make a positive number
                     btrack(nn+tm,13)=latsum(tm)/ngpsum(tm) ! Get mean latitude
                  ENDDO
                  DEALLOCATE(sinsum,cossum,latsum,ngpsum)

                  nn=nn+INT(bstat(label,2))
               ENDIF ! logging
               ! Persistence criterium
               IF ((MAXVAL(ttva(1:ngp))-MINVAL(ttva(1:ngp))+1).LT.persistence) THEN
                  !print*,"Delete label:",label,"from ",MINVAL(ttv(1:ngp)), " to ", MAXVAL(ttv(1:ngp))  
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
      WHERE(arr3(:,:,t2).GT.1) arr3(:,:,t2)=1
      arr2d=arr3(:,:,t2)
      istat = NF90_PUT_VAR(ncido,DatVarIDo,arr2d,start = (/1,1,t1/),count=(/nLong,nLats,1/))
   ENDIF ! IF(tt.GT.(tlen-10)) THEN


ENDDO ! DO tt
istat = NF90_CLOSE(ncid)
istat = NF90_CLOSE(ncido)

IF (logging) THEN 
   OPEN(21,FILE="blockstat_all.txt",ACTION="write",STATUS="unknown",POSITION="rewind")
   OPEN(22,FILE="blockstat.txt",ACTION="write",STATUS="unknown",POSITION="rewind")
   OPEN(23,FILE="blocktracks.txt",ACTION="write",STATUS="unknown",POSITION="rewind")
   OPEN(24,FILE="blocktracks_all.txt",ACTION="write",STATUS="unknown",POSITION="rewind")
   900 FORMAT (2(A8),A5,3(A3),2(A8),A10,  A8,2(A5),2(A9))
   901 FORMAT (2(I8),I5,3(I3),2(I8),F10.3,I8,2(I5),2(F9.3))
   902 FORMAT (A8,2(A7),A10,  A5,3(A3),2(A9),  A11,2(A9))
   903 FORMAT (I8,2(I7),F10.3,I5,3(I3),2(F9.3),I11,2(F9.3))
   WRITE(21,900) 'ID',"Length","YYYY","MM","DD","HH","Start","End","Max_dz","Timx","Lonx","Latx","Lon","Lat"
   WRITE(22,900) 'ID',"Length","YYYY","MM","DD","HH","Start","End","Max_dz","Timx","Lonx","Latx","Lon","Lat"
   WRITE(23,902) 'ID',"Length","Tmstp","Max_dz","YYYY","MM","DD","HH","Lonx","Latx","Area","MeanLon","MeanLat"
   WRITE(24,902) 'ID',"Length","Tmstp","Max_dz","YYYY","MM","DD","HH","Lonx","Latx","Area","MeanLon","MeanLat"
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

SUBROUTINE spoverlapnew(arrover,nLong,nLats,latitude,dx,dy,overlap)

IMPLICIT NONE

   INTEGER,INTENT(IN)                                   :: nLong,nLats
   REAL,INTENT(IN)                                      :: dx,dy,overlap
   REAL,INTENT(IN),DIMENSION(nLats)                     :: latitude
   INTEGER,INTENT(INOUT),DIMENSION(nLong,nLats,3)       :: arrover
   INTEGER                                              :: nn,ii,jj
   REAL(8)                                              :: pi=4.0*ATAN(1.D0)
   REAL(8)                                              :: a1,ax,afrc

   DO nn=2,INT(MAXVAL(arrover(:,:,2))) ! For each contour
      a1=0 ; ax=0 ; afrc=0 
      DO ii=1,nLong
         DO jj=1,nLats
            ! Get area of contour at time tt and tt+1 --> Overlap between both timesteps (forward)
            IF (arrover(ii,jj,2).EQ.nn) a1=a1+ABS(12321.*dx*dy*cos(pi/180*latitude(jj)))
            IF ((arrover(ii,jj,2).EQ.nn).AND.(arrover(ii,jj,3).GE.1)) ax=ax+ABS(12321.*dx*dy*cos(pi/180*latitude(jj)))
            afrc=ax/a1
         ENDDO ! jj
      ENDDO !! ii

      IF (afrc.LT.overlap) THEN
         WHERE(arrover(:,:,2).EQ.nn) arrover(:,:,2)=0
      ENDIF ! IF overlap not big enough
   ENDDO ! nn
END SUBROUTINE spoverlapnew


