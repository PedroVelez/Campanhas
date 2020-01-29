clear all; close all; clc

filedata='Se1507';load(filedata);
load ../DatosCampanha


temmax=24;temmin=2;
salmax=37.025;salmin=34.75;
dens=sw_dens(salts,temps,press)-1000;

%parametros del ts
ts_axis=[salmin salmax temmin temmax];
ts_sigma=[25.:0.5:48.5];

%trozos a cortar
schunks=[salmin-0.2:0.025:salmax+0.2];
tchunks=[temmin-1:0.25:temmax+1];

TS_volumetric(schunks, tchunks,salts,temps,ts_axis,ts_sigma);
title(campanhacode)