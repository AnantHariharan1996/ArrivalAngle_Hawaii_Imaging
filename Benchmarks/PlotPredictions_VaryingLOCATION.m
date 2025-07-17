%% Example Benchmark; plot AA Anomalies as a function of varying TAU
clear; clc; close all;

addpath(genpath('UsefulFunctions'))
load coastlines
FigFolder =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/Stored_ModelSpacePredictions/SummaryFigures/';
RawDatFolder =     '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/Raw_ArrivalAngleDeviations/';
EVID = '11124';
Period = 50;
EVINFO = load([RawDatFolder num2str(Period) 's/' EVID num2str(Period) 's_elon_elat_lon_lat_phidev_phigc' ],'-ascii')
[phvel]=prem_dispersion(1/Period);
cglb = phvel;
spacing = 0.25;
Fname= ['Stored_ModelSpacePredictions/' num2str(Period) 's/ModelSpaceSearch_Store_EVID' ...
    num2str(EVID) '.mat'];

load(Fname)
% get unique pairs of coordinates
zzz(:,1) = Output_Store.Lonstore;
zzz(:,2) = Output_Store.Latstore;
[uniquerows,uniquedx]= unique(zzz,'rows');
UniqueLons = uniquerows(:,1);
UniqueLats = uniquerows(:,2);
UniqueWidths = unique(Output_Store.Widthstore);
UniqueTaus = unique(Output_Store.Taustore);

% Hold Fixed
Fix_Width = UniqueWidths(2);
Scatterer_Lon = -158; Scatterer_Lat = 20;
Scatterer_LON_LIST = linspace(-162,-154,9);
nrows = sqrt(length(Scatterer_LON_LIST));
ncols=nrows;
junk2  = figure(50);
ax = subplot_custom_make(junk2,nrows,ncols,[0.11],[0.25],[0.1 0.9],[0.1 0.9]);
xlist = [-164:0.25:-152];
ylist = [15:0.25:27];
[Ref_XGrid,Ref_YGrid] = meshgrid(xlist,ylist);
EvLon = EVINFO(1,1);
EvLat = EVINFO(1,2);
Scatterer_Tau = 10;

for LonNum = 1:length(Scatterer_LON_LIST)


% Generate Predictions
            scattererlon = Scatterer_LON_LIST(LonNum);
            scattererlat = Scatterer_Lat;
            current_timelag = Scatterer_Tau
            current_width  = Fix_Width;
            % Running forward model! %
        [delta,xgrid,ygrid,ttime_field_perturbed,tau,ttime_field_noscatter,angle,angle_nodiff] ...
            = Get_Arrival_Angle_Residual_GaussianBeam(Period, EvLat,...
            EvLon, scattererlat,scattererlon,current_timelag,current_width,...
            Ref_YGrid(:),Ref_XGrid(:),cglb,spacing);
        
scatter(ax(LonNum),xgrid(:),ygrid(:),20,delta(:),'filled');
title(ax(LonNum),['Diffractor Long. = ' num2str(scattererlon) ' ^\circ'],'fontsize',14)


hold(ax(LonNum),'on')
plot(ax(LonNum),coastlon,coastlat,'linewidth',2,'color','k')
ax(LonNum).XLim = [min(xlist) max(xlist) ];
ax(LonNum).YLim = [min(ylist) max(ylist) ];
clim(ax(LonNum),[-5 5])
ax(LonNum).Box = 'on';
set(ax(LonNum),'linewidth',2)
ax(LonNum).Layer = 'top';
viscircles(ax(LonNum),[scattererlon scattererlat],km2deg(Fix_Width),'linewidth',4)
ax(LonNum).YTickLabels = [];
ax(LonNum).XTickLabels = [];
pbaspect(ax(LonNum),[1 1 1])

IntermediatePt = track2(EvLat,EvLon,scattererlat,scattererlon,100);
IntermediatePts2Use = IntermediatePt(70,:);
quiver(ax(LonNum),IntermediatePts2Use(2),IntermediatePts2Use(1),scattererlon-(IntermediatePts2Use(2)),scattererlat-IntermediatePts2Use(1),'linewidth',2,'MaxHeadSize',0.4,'color','magenta')


end
set(gcf,'Position',[1865 158 654 621])
barbar=colorbar(ax(LonNum));
barbar.Location = 'southoutside';
barbar.Position = [ 0.125  0.07   0.745 0.0164];
barbar.FontSize = 14;
ylabel(barbar,'Arrival Angle Deviation (degrees)','FontWeight','bold')
colormap(turbo)
sgtitle('Impact of Varying Diffractor Location on Arrival Angles','fontsize',20)
saveas(gcf,'Benchmarks/VaryingLocExample.png')