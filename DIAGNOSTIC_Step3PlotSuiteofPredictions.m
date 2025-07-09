% Interrogate Model Space

% Make figure showing AA predictions for different model parameters,
% for a different event
clear; clc; close all; 
addpath(genpath('UsefulFunctions'))
load coastlines
FigFolder =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/Stored_ModelSpacePredictions/SummaryFigures/';

EVID = '11124';
Period = 50;
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
Fix_Tau = 10;

%%% 
junk2  = figure(50);
ax = subplot_custom_make(junk2,5,5,[0.11],[0.25],[0.1 0.9],[0.1 0.9]);

for modelnum =1:length(ax(:))

    current_lon = randsample(UniqueLons,1);
    current_lat = randsample(UniqueLats,1);
    pred_dx = find(Output_Store.Lonstore == current_lon ...
        & Output_Store.Latstore == current_lat ...
        & Output_Store.Widthstore == Fix_Width ...
       & Output_Store.Taustore == Fix_Tau)

   scatter(ax(modelnum),Output_Store.xgrid,Output_Store.ygrid,10,Output_Store.ModelSpaceSearch_Store(pred_dx,:),'filled')
hold(ax(modelnum),'on')
plot(ax(modelnum),coastlon,coastlat,'linewidth',2,'color','k')
ax(modelnum).XLim = [min(Output_Store.xgrid) max(Output_Store.xgrid) ]
ax(modelnum).YLim = [min(Output_Store.ygrid) max(Output_Store.ygrid) ]
clim(ax(modelnum),[-5 5])
ax(modelnum).Box = 'on'
set(ax(modelnum),'linewidth',2)
ax(modelnum).Layer = 'top'
viscircles(ax(modelnum),[current_lon current_lat],km2deg(Fix_Width))
ax(modelnum).YTickLabels = [];
ax(modelnum).XTickLabels = [];
pbaspect(ax(modelnum),[(max(Output_Store.xgrid)-min(Output_Store.xgrid))/(max(Output_Store.ygrid) - min(Output_Store.ygrid) ) 1 1])


end
colormap(turbo)

sgtitle({['Time Delay \tau = ' num2str(Fix_Tau) 's,',' Width = ' num2str(Fix_Width) ' km, Period = ' num2str(Period) 's'], 'Color Ranges from -5 deg to + 5 deg Deviation from Great-Circle Path'},'fontsize',26,'fontweight','bold')
set(gcf,'position',[1571 -41 1442 906])
saveas(gcf,[FigFolder 'EVID' num2str(EVID) '_FixedTau_' '.png'])

