%% More Benchmarking. 
clear; clc; close all;
% generate arrival angles on synthetic traveltime data,
% and compare it to great-circle calculations. 

Evlo=0; Evla=0; phvel=4;
stalons = [-180:1:180];
stalats = [-90:1:90];
[LONGRD,LATGRD] = meshgrid(stalons,stalats);
[dist2stns,az2stns] = distance(LATGRD,LONGRD,Evla,Evlo);
dist2stns_km = deg2km(dist2stns);
ttime = dist2stns_km./phvel;
[fx,fy]=gradient(ttime,deg2km(1),deg2km(1)); %dt/dtheta and dt/dphi
sincolat = sind(90-LATGRD);
arrival_angle = rad2deg( atan2((fx./sincolat),fy)    );
%arrival_angle(arrival_angle<0)=arrival_angle(arrival_angle<0)+360;
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
scatter(LONGRD(:),LATGRD(:),5,arrival_angle(:),'filled')
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
scatter(LONGRD(:),LATGRD(:),5,az2stns(:)-arrival_angle(:)-180,'filled')
hold on
grid on; box on;
plot(coastlon,coastlat,'linewidth',2,'color','k')
barbar=colorbar;
ylabel(barbar,'Source-Station Azimuth  (deg)')
title('Difference Between Calculation and Great circle Prediction and minus 180')
set(gca,'fontsize',18,'fontweight','bold')
caxis([-10 10])

set(gcf,'position',[84 146 1382 684])
saveas(gcf,'BenchmarkOutput.png')