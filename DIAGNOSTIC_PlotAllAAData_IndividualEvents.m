%% Diagnostic step; plot misfit surfaces for each event. 

clear; clc; close all; load coastlines; 
HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';
PredictionsDir =  [HomeDir 'Stored_ModelSpacePredictions/'];
ObservationsDir = [HomeDir 'Raw_ArrivalAngleDeviations/'];
SummaryMisfitDir = [PredictionsDir 'SummaryMisfitStore/'];
FigFolder =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/Stored_ModelSpacePredictions/SummaryFigures/';
Cmap2Use =     '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/UsefulFunctions/roma.cpt';
Periodlist = [50 66.6667 80 100];

for  Period=Periodlist
    ModelStorageFolder = [PredictionsDir num2str(Period) 's/'];
    ObservationsDir_CurrentPeriod = [ObservationsDir num2str(Period) 's/'];
    MisfitSurfaceFigs_CurrentPeriod = [PredictionsDir num2str(Period) 's_Figs/'];


Step2b_SetupParameters

% For EVERY event, plot the observations' implied geographical weights


ID_Fname = [PredictionsDir 'IDLIST_' num2str(Period) 's.txt'];
IDList = load(ID_Fname)
numevts = length(IDList);

junk2 = figure(149)
nrows = ceil(sqrt(numevts));
ncols =nrows;
ax = subplot_custom_make(junk2,nrows,ncols,[0.05],[0.05],[0.05 0.95],[0.05 0.95]);
for eventnum = 1:length(IDList)
    EVID = num2str(IDList(eventnum));
    % Load the observations.    
    RawObs_Fname = [ObservationsDir_CurrentPeriod EVID num2str(Period) 's_elon_elat_lon_lat_phidev_phigc'];
    AngObsInfo =  load(RawObs_Fname,'-ascii');
    Elon = AngObsInfo(1,1); Elat = AngObsInfo(1,2);
    Obs_Xgrid = AngObsInfo(:,3); Obs_Ygrid = AngObsInfo(:,4); AngleResid = AngObsInfo(:,5);
    Obs_Weights = AngObsInfo(:,8);
     RawStationData_Fname = [ObservationsDir_CurrentPeriod EVID num2str(Period) 's_slon_slat_stt'];
    RawStnInfo =  load(RawStationData_Fname,'-ascii');  
    RawSlon = RawStnInfo(:,1);     RawSlat = RawStnInfo(:,2);     RawStt = RawStnInfo(:,3);


scatter(ax(eventnum),RawSlon,RawSlat,100,[0.5 0.5 0.5],'^','filled')

hold(ax(eventnum),"on")
scatter(ax(eventnum),Obs_Xgrid,Obs_Ygrid,20,AngleResid,'filled')
scatter(ax(eventnum),RawSlon,RawSlat,100,[0.5 0.5 0.5],'^','filled','linewidth',2)

plot(ax(eventnum),coastlon,coastlat,'color','k','LineWidth',2)
text(ax(eventnum),-158, 23,['Event ' EVID],'fontsize',16,'FontWeight','bold')
   
% Quote and display QC Parameters
% text(ax(eventnum),-162, 23,['RMS = ' num2str(round(current_RMS,2)) 's'],'fontsize',16)
%   text(ax(eventnum),-162, 20,['N = ' num2str(round(length(Obs_Xgrid),2)) ],'fontsize',16)

%box on
set(ax(eventnum),'linewidth',2,'box','on','XTickLabel',[],'YTickLabel',[])
   ax(eventnum).XLim=[-159.5 -152.5];
     ax(eventnum).YLim=[17 25];
ax(eventnum).CLim = [-5 5];
colormap(ax(eventnum),jet(20))
 cptcmap(Cmap2Use,ax(eventnum),'ncol',20);

end

%%% 


set(gcf,'position',[1637 -77 1134 956])
currnum = eventnum;

for eventnum = currnum+1:1:nrows*ncols

ax(eventnum).Color = 'none'
ax(eventnum).Box = 'off'
ax(eventnum).XTick = [];
ax(eventnum).YTick = [];
ax(eventnum).XColor = 'none'
ax(eventnum).YColor = 'none'
ax(eventnum).FontSize = 18

if eventnum == currnum+1
% add colorbar
barbar = colorbar(ax(eventnum),'location','west')
 cptcmap(Cmap2Use,ax(eventnum),'ncol',20);

ax(eventnum).CLim = [-5 5];
ylabel(barbar,{'Arrival Angle','Deviation (˚)'}, 'fontsize',18)
end
end


text(ax(eventnum),0.5,0.5,[num2str(Period) ' s'],'fontsize',24,'fontweight','bold')

saveas(junk2,[FigFolder num2str(Period) '_sAllObs.png'])
end

