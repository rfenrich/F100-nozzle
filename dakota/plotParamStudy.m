clc; clear all; close all;

A = importdata('TurbofanParamStudy.dat');
A = A.data;

alt = A(:,1); mach = A(:,2); sfc = A(:,3);

[ALT MACH] = meshgrid(sort(unique(alt)), sort(unique(mach)));
SFC = griddata(alt, mach, sfc, ALT, MACH);

figure;
contourf(ALT, MACH, SFC,30);
shading flat;
xlabel('Altitude, [ft]');
ylabel('Mach Number');
h = colorbar;
ylabel(h,'sfc, [lb/lbf/hr]');
