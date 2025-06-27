%% Compare speed of distance functions

clear; clc; close all; 

tic  
[ KMLEN,AZ] = distance_deg_nomap(0,0,80,80)
toc

tic  
[ KMLEN,AZ] = distance(0,0,80,80)
toc

tic 
[tracklat,tracklon] = track1(0,0,80,180,[],[],1000);
toc

tic
toc