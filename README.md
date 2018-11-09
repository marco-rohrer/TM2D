# TM2D - Tibaldi Molteni 2-dimensional blocks
2D dimensional blocking algorithm (see Rohrer et al., 2018; Rohrer et al., submitted to Tellus A; Rohrer et al., in prep.). Computed blocking catalogues for a multitude of dataset are available on request (marco.rohrer@giub.unibe.ch ; mrohrer87@gmail.com).

The algorithm detects potentially blocked regions in a gridded dataset. As of now three different methods are implemented: TM2D, Z500anom, VAPV. The name of this repository is thus a little out of date, as more methods were added later. After potentially blocked areas are detected, the algorithm applies a persistence criteria and an overlap criteria that ensures that a block does not move too fast between two time steps. In its standard setting the persistence criterion is set 20 time steps, which corresponds to 5 days for 6-hourly data. The standard setting for the overlap criterion is 0.7 or 70% overlap between two time steps (i.e. A<sub>t</sub> ⋂ A<sub>t+1</sub> ≥0.7*A<sub>t</sub>) . This is almost identical to the tracking approach in Schwierz et al. (2004). 
On a more technical note, the algorithm was programmed in a way to detect blockings also on very long datasets such as CCC400 with 400 years of data. One large NetCDF file is created that is loaded sequentially into the blocking algorithm. This approaches enables to run the software with a relatively low RAM footprint, although very large input files (e.g. a complete CCC400 member is 40 GB large). Instead of loading the whole dataset at once, only approximately one year is loaded into the algorithm, reducing the memory footprint by a factor of roughly 400. As soon as a time step is not used anymore, as blockings can only be detected forward in time but not backward, the next time step is loaded until the end of the file is reached. This approach ensures that every detected blocking has an unique ID and evades problems arising if the datasets is examined in e.g. yearly chunks (i.e. time overlap needed to detect blockings at the beginning and end of the chunk). 

## Supported methods
Here, the three supported methods are briefly described. For more information refer to the Rohrer et al. papers. 

### TM2D
TM2D combines the approach of Scherrer et al. (2006), who uses a two-dimensional extension of the original Tibaldi and Molteni (1990) algorithm based on the detection of reversals of the meridional 500 hPa geopotential heigt gradient. Following two criteria are required to flag a grid point as potentially blocked:
a) GPH gradient towards the pole:       GPHG<sub>P</sub>=  (Z500<sub>φ+dφ</sub> - Z500_<sub>φ</sub> )/dφ < -10 gpm/(°lat)
b) GPH gradient towards the equator:	GPHG<sub>E</sub> =  (Z500<sub>φ</sub>) - Z500<sub>φ-dφ</sub>)/dφ  > 0  gpm/(°lat)
The difference between the central and the polarward or equatorward grid point dφ is given by INT((15°)/dy)*dy (where dy the horizontal resolution in latitudinal direction. Therefore, if a 2° resolution is chosen, this results in  dy=INT((15°)/(2°))*2°= 14° and the latitude φ varies from 36° to 76° in 2° intervals. 


### Z500anom
This approach detects potential blocks as regions with a positive Z500 anomaly (Z500*), similar to Dole and Gordon (1983). Here, we use regions with a Z500 anomaly larger than 200 m (again derived by subtracting the 30 day running mean climatology). Please note that the algorithm does NOT compute the anomaly on its own, so the input file needs to already contain Z500 anomalies. Of course also other variables than Z500 can potentially be used.

#### VAPV
Similar to the definition in Schwierz et al. (2004), Pfahl & Wernli (2012b) and Lenggenhager et al. (in revision) blocks are derived from the 150–500 hPa vertically averaged PV (VAPV). First the VAPV anomaly is computed by subtracting the 30 day running mean climatology of the corresponding time step of the year. To smooth out high-frequency variability an additional two day running mean filter is applied before detecting regions with PV <=-1.3 in the Northern Hemisphere and >=1.3 in the Southern Hemisphere. Note that internally values in the Southern Hemisphere are multiplied by (*-1) and the blocking statistics files are not corrected for that! I.e. a block with an anomaly of -2 PVU in the Southern Hemisphere has an anomaly of +2 PVU in reality.

## Pre- and postprocessing
This folder contains 3 sample bash scripts that take care of the pre- and postprocessing. Preprocessing includes remapping and calculating anomalies. Note that for VAPV and Z500anom the exact choice of the anomaly (i.e. compare it to a fixed 30 year climatology vs. to a transient 30-year running mean climatology) may change the results! The algorithm expects one large file, so at the end we merge the whole input dataset to one large NetCDF file. The file should have an absolute time axis is required, otherwise the timestamp is not correctly computed. An absolute time axis can be created by the command: cdo **-a** copy trel.nc tabs.nc
Postprocessing includes moving files to a predefined folder and computing seasonal, yearly block climatologies. Do not forget to change the representation type from INTEGER (as we only have labelled integers in the output NetCDF file) to FLOAT to represent fractions (i.e. cdo **-b F32** timmean blocks.nc blockfreq.nc). Additionally text files with blocktracks and block statistics are created. blocktracks.txt and blockstat.txt contain data of all valid blocks, blocktracks_all.txt and blockstat_all.txt contain information about discarted blocks that are NOT included in the output NetCDF. 

## Command line interface
TM2D uses a command line interface that can change the most important parameters without the need to actually change and recompile any code. Here a list of all options is given:  
--logging=BOOL : Use this option to deactivate the log files that saves statistics about the blocks [T]  
--infile=CHAR: REQUIRED, use this parameter to indicate the file that is read in, e.g. /home/user/Z500.nc  
--invar=CHAR : REQUIRED, use this parameter to indicate which variable in infile is read on, e.g. Z500  
--outfile=CHAR : REQUIRED, use this parameter to indicate the file to which the blocks are written, e.g. /home/user/blocks.nc  
 --morehelp [Print even more help]  
 --infile=CHAR [Indicate input file name - REQUIRED]  
 --invar=CHAR [Indicate required variable name of infile - REQUIRED]  
 --outfile=CHAR [Indicate output file name - REQUIRED]  
 --mode=CHAR [Indicate blocking mode: {TM2D,Z500anom,VAPV}]  
 --persistence=INT [# timesteps for persistence, required for TM2D,Z500anom,VAPV][20]  
 --overlap=FLOAT [required spatial overlap, required for TM2D,Z500anom,VAPV][0.7]  
 --dlat=INT [degrees between central and N/S grid point for TM2D [15]  
 --dzp=INT [required gph gradient towards pole [m/°lat]for TM2D [-10]  
 --dze=INT [required gph gradient towards equator [m/°lat] for TM2D [0]  
 --minlat=INT [minimum latitude towards equator for TM2D [35]  
 --anom=INT [minimum geopotential height anomaly [m] for Z500anom [200]  
 --vapvmin=FLOAT [minimum VAPV treshold [PVU] for VAPV [-1.3]  
 -nologging [deactivate logging of block characteristics, increases performance]  
 -nooverlap [deactivate overlap criterium]  
 -nopersistence [deactivate persistence criterium]  
 -nc3 [creates an uncompressed NetCDF3 file instead of a zipped NetCDF4 file  
 -labels [activate ascencding numbering of blocks instead of binary field]  

### Installation
You need a working NetCDF Fortran library installed on your computer with HDF5 support, otherwise use the -nc3 flag to omit the creation of a compressed NetCDF-4 file. Additionally, gfortran (any recent version, > 4.3? should do) must be installed. 
The tm2d.make script compiles the tm2d fortran code and should automatically link the NetCDF-Fortran library that is found first in our $PATH variable. 
The now created executable can be called like a normal program, i.e. ./tm2d

### Problems and Bugs
In case of bugs, either send me a bug report with as much information as possible (mrohrer87@gmail.com) or fix it yourself :) In that case let me know, what failed and how you fixed it. 

### Example
./tm2d --infile=/home/marco/Z500.nc --invar=Z --outfile=/home/marco/blocks.nc --mode=TM2D

This produces TM2D blocks using the file /home/marco/Z500.nc as its input and extracts variable Z from this file as its input variable (you can check the variable name with ncdump, in case you don't know). After the block detection a NetCDF file blocks.nc is created. Additionally the blockstat and blocktrack text files are generated.

./tm2d --infile/home/marco/VAPV.nc --invar=PV --outfile=/home/marco/PVblocks.nc --mode=VAPV --vapvmin=-1.0 -labels -nologging -nc3

This produces VAPV blocks using the variable PV from file /home/marco/VAPV.nc. Some additional parameters are given that alter the behaviour of the algorithm. The --vapvmin=-1.0 flag changes the standard setting of -1.3 PVU as the threshold to detect a potential block to -1.0, i.e. weakens the PV criterium. The -lables flag indicates that the outputted NetCDF file contains not a binary [0/1] field of nonblocked/blocked gridpoints but that the labels (i.e. the ID in the blocktracks and blockstat files) are retained in the file. -nologging disabels the creation of the blocktracks and blockstat files, this can increase performance considerably for large datasets but only creates the output NetCDF file. Finally the -nc3 option disables the NetCDF-4 file compression and creates a classic NetCDF-3 file. This file will be massively larger, so only use that option in case you run into problems with NetCDF-4/HDF5. 

### Remarks
Note that this version is not thoroughly checked yet, so bugs can still occur. Without extra options there should be no problem, but not every combination of extra options was checked for mode of calculation.

### Possible future extensions of the algorithm
- Deactivate switch for blockstat_all.txt and blocktracks_all.txt
- Better documentation of code
- Switch for tlen (as a --maxpersistence switch) to indicate the maximum length a block can possibly have
- Insert a criteria for the size of blocks (like in the original Schwierz algorithm)
- Code clean up
- Arrange the overlap and contouring subroutines as external modules as these algorithms may be handy for other programs too


Marco Rohrer (marco.rohrer@giub.unibe.ch ; mrohrer87@gmail.com), Nov 2018

