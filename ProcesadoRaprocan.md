# CTD and LADCP data procesing

This manual describes how to process the SBE911+ and LADCP data during the Raprocan cruises.

The data is organized in the [Campanha] folders as follows:
```
[Campanha]/Analisis
     --- /Hidrografia
[Campanha]/CTD
     --- /Raw {.hex data}
     --- /Mat {.mat data}
     --- /Cnv {.cnv data}       
[Campanha]/LADCP
     --- /Raw {.cnv data}       
     --- /Procesados {.cnv data}                        
     --- /Procesados/profiles {processed data}                                         
     --- /Procesados/plots {processed data}
```

## Hidrografia
The CTD data processing is carried out in the [Campanha]/Analisis/Hidrografia folder

First, yo have to  set up the configuration file for the cruise:

* Edit *DatosCampanha.m* in [Campanha], and run it from the Analisis folder.
```
campanha='Raprocan1810'; % Cruise name
campanhacode='Ra1810'; % Short cruise name
dirdata=''[...]/[Campanha]/CTD';
Also the geographical area, coast file,...
```

### Convert cnv files to -Mat
* Copy the .cnv to process into the [Campanha]/CTD/Cnv into the [Campanha]/Analisis/Hidrografia folder

Call *ProcesaSbe911.m*
This script create de .mat files with the variables. and move it to [...]/[Campanha]/CTD/Mat' . Option (2) will  move the .mat file to [Campanha]/CTD/Mat,but it will not delete the .cnv file

Set Automatico=1 to process all stations at once. If set to 0, stations will be processed one by one, allowing to visualize the profiles. Once the process has been done, the .cnv files can be deleted from the Hidrografia folder

### Create a mat file with the selected cast files
Call CreaSecCTD.m (from Hidrografia folder)  and indicate the cast files you want in a section.
There are 2 options: create a unique .mat files with all the sections (for instance [1:1:24 901:1:910]) or create one .mat file per sections. The second option is better for plotting variables sections (with secciones.m). Usually the name is Ra1903.mat (all CTD casts) or Ra1903Norte.mat for the 1:1:24 CTD casts.

### EvolucionDiff.m
Monitor the time evolution of the difference between the dual T and C (for salinity) sensors. To be called from Hidrografia folder. Edit the datafile indicating the name of .mat file with the data.

### TS diagram for the selected .mat file
Call *DiagramaTS.m* from Hidrografia folder.
Edit the name of the datafile (.mat file for each section or unique .mat file) and edit zones limits and names (if sections was changed) . The .mat used is the one created CreaSecCTD.m

### Vertical sections
Call *Secciones.m* from Hidrografía folder.
Edit line where the data is loaded:
```
Data=load('Ra1903Norte');
```
For a new section, we need to create a bathymetric profile first, calling CreaBatimetriaRadial.m from Hidrografia folder
Edit also the opciones ('z' for zonal and 'm' for meridional), Data, SVOpciones.tipo and SVOpciones.batymetry (.m bathymetry file just created with CreaBatimetriaRadial.m)

### Compare with previous CTD stations.
ComparaEstacionesPrevias.m

## Transport estimates.
Transport estimates are computed in [Campanha]/Analisis/Transporte

First, it is necessary to modify the configuration files
* Rename and modify the [..]/[Campanha]/Analisis/Transporte/sig_[campanha]_gamma.dat
* Modify the cntrl_[SectionName].dat file for each section where the transport would be computed.
* Rename, get_stas_[campanha].m and Modify flname of the new CTD files.
* Rename, mk_circ2_[campanha].m and Modify Call to get_stas_[campanha]

### Compute transports
* Rename Circ3s_[campanha].m and modify:
```
output_dir='[...]/[campanha]/Analisis/Transporte/';
layer_file='sig_[campanha]_gamma.dat';
Rename calls to mk_circ2_[campanha].m
```
* call *Circ3s_[campanha].m*

### Visualize Layer transport
* Edit g_transporte_capa_[Section].
Change transport file (line 6)
Comment lines 9 to 11 if reference level is not calculated yet. If it's the case, uncomment line 12.

### Visualize accumulated transport
* Edit g_transporte_acumulado_[Section].
```
stations=[11:24]; %Stations in the section
iCapaSuperior=1:3;
iCapaIntermedia=4:6;
iCapaProfunda=7:12;
Edit the call to get_stas_[campanha]
```
Comment lines 18 to 20 if reference level is not calculated yet. If it's the case, uncomment line 21.

## LADCP data proccesing

Copy folder [Campanha] in [...]/ProcesadoLADCP/ProcesadoLADCP10/[campanha]

Edit the startup.m file and set up the path for the LADCP10 folder: [...]/ProcesadoLADCP/ProcesadoLADCP10/[campanha]

Edit cruise name in [...]/ProcesadoLADCP/ProcesadoLADCP10/[campanha]/m/
* cruise_params.m
* prepctdprof.m
* prepctdtime.m (line 73,74)[Cruise name and year]
* prepnav.m (line 50)

### Process CTD data for ctdprof and ctdtime
* Copy .hex files to /Users/pvb/Intercambio/Raprocan/Cnv
* Edit /Users/pvb/Intercambio/Raprocan/raw_ctdprof/batch_ctdprof.txt
    change the folder, CTD configuration file and cruise name
* Edit /Users/pvb/Intercambio/Raprocan/raw_ctdprof/DatCnv.psa
    change the configuration file
* Edit /Users/pvb/Intercambio/Raprocan/raw_ctdtime/batch_ctdtime.txt
    change the folder, CTD configuration file and cruise name
* Edit /Users/pvb/Intercambio/Raprocan/raw_ctdtime/DatCnv.psa
    change the configuration file
* Open Virtual box and run scripts in /Users/pvb/Intercambio/Raprocan/sbe911batchLeeME.txt and run the SEeaBird data Procesing

* Copy from .../raw_ctdtime to /ProcesadoLADCP10/[campanha]/data/raw_nav and ctdtime
* Copy from .../raw_ctdprof to /ProcesadoLADCP10/[campanha]/data/raw_ctdprof

### Process LADCP datafile
* Copy raw LADCP data into /ProcesadoLADCP10/[campanha]/data/raw_ladcp
each stations should be in a folder with the following name:
```
Master /[Campanha]/data/raw_ladcp/NNN/NNNdn000.000  [NNN is the station number]
Slave /[Campanha]/data/raw_ladcp/NNN/NNNup000.000  [NNN is the station number]
```

*open matlab, and go to the ...]/ProcesadoLADCP/ProcesadoLADCP10/[campanha] folder
```
process_cast(NumeroEstacion)
```
## Satelite data
Bajar datos de satelite del

VIIRSL3m
%https://oceandata.sci.gsfc.nasa.gov/VIIRS-SNPP/Mapped/8-Day/4km/sst/
%https://oceandata.sci.gsfc.nasa.gov/VIIRS-SNPP/Mapped/8-Day/4km/chlor_a/

y el AquaL3m
https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/Mapped/8-Day/4km/chlor_a/
https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/Mapped/8-Day/4km/sst/

### SST and Dynamic topography
Madt -> Mapa de MADT
SST > Mapas de SST
