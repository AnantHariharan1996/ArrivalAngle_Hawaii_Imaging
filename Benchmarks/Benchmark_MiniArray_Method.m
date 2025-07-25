%% Solve for Arrival Angles 
%% Using a Mini-array approach
fundir =     '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/UsefulFunctions';
addpath(genpath(fundir))
clear; clc; close all;
load coastlines
xlist= [-120:0.8:-80];
ylist= [28:0.8:48];
Phvel = 4;
%% Inputs
info = readtable('STAFILE_USARRAY.csv');
Sta_Lons = info.Var3;
Sta_Lats=info.Var2;
EvLon = [-130];
EvLat = [10];
distlist=distance(EvLat,EvLon,Sta_Lats,Sta_Lons);
distlist_km=deg2km(distlist);
L_Tol = 1;
TTime = distlist_km./Phvel;
Min_N_Stations=4;


%%

[ArrivalAngleList,Best_LocalPhVelList] = ...
    GetArrivalAngles_Event_MiniArray(Sta_Lons,Sta_Lats,TTime,L_Tol,EvLat,EvLon,Min_N_Stations);
size
[JUNK,STRAIGHTANGLE] = distance(Sta_Lats,Sta_Lons,EvLat,EvLon);
[ fx,fy,angle,xgrid,ygrid,tgrid2 ] = Get_arrival_angle( EvLat,EvLon,Sta_Lats,Sta_Lons,TTime,0.5 );
figure(1)
subplot(2,2,1)
plot(coastlon,coastlat,'linewidth',2,'color','k')
hold on
scatter(EvLon,EvLat,200,[1 0 0],'pentagram','filled')
scatter(Sta_Lons,Sta_Lats,20,[0 0 1],'filled','^')
xlim([-140 -67.5])
ylim([10 50])
title('Source-Receiver Geometry')
set(gca,'fontsize',18)

subplot(2,2,4)
scatter(Sta_Lons,Sta_Lats,20,ArrivalAngleList,'filled')
barbar=colorbar;
caxis([10 80])
title('Arrival Angle from Mini-Array Approach')
set(gca,'fontsize',18)
xlim([-130 -67.5])
ylim([25 50])
ylabel(barbar,'Arrival Angle (degrees)','fontsize',18)

subplot(2,2,3)
scatter(xgrid,ygrid,20,angle,'filled')
barbar=colorbar;
caxis([10 80])
title('Arrival Angle from Traveltime Gradient')
set(gca,'fontsize',18)
xlim([-130 -67.5])
ylim([25 50])
ylabel(barbar,'Arrival Angle (degrees)','fontsize',18)

subplot(2,2,2)
scatter(Sta_Lons,Sta_Lats,20,Best_LocalPhVelList,'filled')
barbar=colorbar;
ylabel(barbar,'Local Phase Velocity (km/s)','fontsize',18)
caxis([3.975 4.025])
colormap(turbo)
title('Local Phase Velocity from Mini-Array Approach')

xlim([-130 -67.5])
ylim([25 50])
sgtitle('Benchmark of Arrival-Angle Estimation Approaches: c = 4 km/s','fontsize',30,'fontweight','bold')
set(gca,'fontsize',18)
set(gcf,'position',[52 58 1248 807])
saveas(gcf,'Benchmark.jpg')

%%% Get Difference Between Approaches
Interped_GradAngles = griddata(xgrid,ygrid,angle,Sta_Lons,Sta_Lats);

figure()
scatter(Sta_Lons,Sta_Lats,20,Interped_GradAngles'-ArrivalAngleList,'filled')
barbar=colorbar;
caxis([-10 10])
title('Interpred Grad Angle - Beamforming Angle')
set(gca,'fontsize',18)
xlim([-130 -67.5])
ylim([25 50])
ylabel(barbar,'Arrival Angle Dev (degrees)','fontsize',18)
colormap(jet)

