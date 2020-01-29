# CTD and LADCP data procesing

This manual describes how to process the SBE911+ and LADCP data during the Raprocan cruises.

The data is organize in the [Campanha] folders as follows:
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

* Edit *DatosCampanha.m* in [Campanha]
```
campanha='Raprocan1810'; % Cruise name
campanhacode='Ra1810'; % Short cruise name
dirdata=''[...]/[Campanha]/CTD';
Also the geographical area, coast file,...
```

### Convert cnv files to -Mat
* Copy the .cnv to process into the [Campanha]/CTD/Cnv into the [Campanha]/Analisis/Hidrografia folder

* Call *ProcesaSbe911.m*
This script create de .mat files with the variables and move it to [...]/[Campanha]/CTD/Mat'
option (2) will  move the .mat file to [Campanha]/CTD/Mat,but it will not delete the .cnv file

### Create a mat file with the selected cast files
Call CreaSecCTD.m and indicate the cast files you want in a section. Usually the name is Ra1903.mat (all CTD casts) or Ra1903Norte.mat for the 1:1:24 CTD casts.

### EvolucionDiff.m
Monitor the time evolution of the difference between the dual T and C (for salinity) sensors.
Edit the datafile

### TS diagram for the selected .mat file
DiagramaTS.m
Use the .mat file created with CreaSecCTD.m
Edit the datafile

### Vertical sections
* Secciones.m
Edit line where the data is loaded, for each section to be carried out.
```
Data=load('Ra1903Norte');
```

### SST and Dynamic topography
Madt -> Mapa de MADT
SST > Mapas de SST

### Compare with previous CTD stations.
ComparaEstacionesPrevias.m

## Transport estimates.
Transport estimates are computed in [Campanha]/Analisis/Transporte

First, it is necesary to modify the configuration files
* Rename and modify the [..]/[Campanha]/Analisis/Transporte/sig_[campanha]_gamma.dat
* Modify the cntrl_[SectionName].dat file for each section where the transport would be computed.
* Rename, get_stas_[campanha].m and Modify
    flname of the new CTD files.
* Rename, mk_circ2_[campanha].m and Modify
    Call to get_stas_[campanha]

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

### Visualize accumulated transport
* Edit g_transporte_acumulado_[Section].
```
stations=[11:24]; %Stations in the section
iCapaSuperior=1:3;
iCapaIntermedia=4:6;
iCapaProfunda=7:12;
Edit the call to get_stas_[campanha]
```

## LADCP data procesing

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

*open matlab,and go to the ...]/ProcesadoLADCP/ProcesadoLADCP10/[campanha] folder
```
process_cast(NumeroEstacion)
```
