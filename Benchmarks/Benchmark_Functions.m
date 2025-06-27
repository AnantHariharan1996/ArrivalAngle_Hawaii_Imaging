%% More Benchmarking. 
addpath(    '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/')
clear; clc; close all;
% generate arrival angles on synthetic traveltime data,
% and compare it to great-circle calculations. 

Evlo=0; Evla=0; phvel=4;
stalons = [-180:1:180];
stalats = [-90:1:90];
[LONGRD,LATGRD] = meshgrid(stalons,stalats);
[dist2stns,az2stns] = distance(LATGRD,LONGRD,Evla,Evlo);

[dist2stns,az2stns] = distance(LATGRD,LONGRD,Evla,Evlo);
dist2stns_km = deg2km(dist2stns);
ttime = dist2stns_km./phvel;

 [ fx,fy,angle,xgrid,ygrid,tgrid2 ] = Get_arrival_angle( Evla,Evlo,LATGRD(:),LONGRD(:),ttime(:),1 );

angle(angle>180)=angle(angle>180)-360;
load coastlines

%%%% Plotting below here
figure
subplot(2,2,1)
scatter(LONGRD(:),LATGRD(:),5,ttime(:),'filled')
hold on
scatter(Evlo,Evla,150,'pentagram','filled')
hold on
grid on; box on;
plot(coastlon,coastlat,'linewidth',2,'color','k')
barbar=colorbar;
ylabel(barbar,'Traveltime (s)')
title('Traveltime')
set(gca,'fontsize',18,'fontweight','bold')

subplot(2,2,3)
scatter(xgrid(:),ygrid(:),5,angle(:),'filled')
hold on
grid on; box on;
plot(coastlon,coastlat,'linewidth',2,'color','k')
barbar=colorbar;
ylabel(barbar,'Arrival Angle from ttime gradient  (deg)')
title('Arrival Angle from Ttime Gradient (deg)')
set(gca,'fontsize',18,'fontweight','bold')

subplot(2,2,4)
scatter(LONGRD(:),LATGRD(:),5,az2stns(:),'filled')
hold on
grid on; box on;
plot(coastlon,coastlat,'linewidth',2,'color','k')
barbar=colorbar;
ylabel(barbar,'Source-Station Azimuth (deg)')
title('Azimuth from Event to Station using ''distance'' function')
set(gca,'fontsize',18,'fontweight','bold')

subplot(2,2,2)
scatter(LONGRD(:),LATGRD(:),5,az2stns(:)-griddata(xgrid(:),ygrid(:),angle(:),LONGRD(:),LATGRD(:))-180,'filled')
hold on
grid on; box on;
plot(coastlon,coastlat,'linewidth',2,'color','k')
barbar=colorbar;
ylabel(barbar,'Source-Station Azimuth  (deg)')
title('Difference Between Calculation and Great circle Prediction and minus 180')
set(gca,'fontsize',18,'fontweight','bold')
caxis([-10 10])

set(gcf,'position',[84 146 1382 684])