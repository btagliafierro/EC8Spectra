clear; close all; clc


g=9.81;

ec8=designSpectrum;   
ec8.behaviourFactor=1.5;   
ec8.importanceFactor=1.0;
ec8.type    = 'Type 1';
ec8.soil    = 'B';
ec8.ag=0.261*g;

figure(1); hold on; box on
plot(ec8.time,ec8.pseudoAcc/g,'b')
axis([0 4 0 inf])
