%% Step : Plot Misfit Surfaces
% What does this script do?
% Make ONE figure that has the misfit surfaces for ALL events. 


clear; clc; close all; load coastlines; 
HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';
PredictionsDir =  [HomeDir 'Stored_ModelSpacePredictions/'];
ObservationsDir = [HomeDir 'Raw_ArrivalAngleDeviations/'];
SummaryMisfitDir = [PredictionsDir 'SummaryMisfitStore/'];
FigFolder =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/Stored_ModelSpacePredictions/SummaryFigures/';

Periodlist = [66.6667 80 100];

for  Period=Periodlist
    ModelStorageFolder = [PredictionsDir num2str(Period) 's/'];
    ObservationsDir_CurrentPeriod = [ObservationsDir num2str(Period) 's/'];
    MisfitSurfaceFigs_CurrentPeriod = [PredictionsDir num2str(Period) 's_Figs/'];

    fname = [SummaryMisfitDir '/MisfitSurfaces_' num2str(Period) 's.mat'];
    load(fname)

Step2b_SetupParameters

% For EVERY event, plot the observations, the 
% best-fitting model, and best-fitting misfit cross-sections
% as a function of width position and tau. 
StoreAllMisfits_L1 = MisfitSurfaceSummary.StoreAllMisfits_L1;
StoreAllMisfits_L2 = MisfitSurfaceSummary.StoreAllMisfits_L2;
PhVelList = MisfitSurfaceSummary.PhVelList;
RMSList = MisfitSurfaceSummary.RMSList;

 Lonstore = MisfitSurfaceSummary.Lonstore;
Latstore=  MisfitSurfaceSummary.Latstore;
Widthstore= MisfitSurfaceSummary.Widthstore;
 Taustore= MisfitSurfaceSummary.Taustore;
EVIDLIST=    MisfitSurfaceSummary.EVIDLIST;
VR_Store = MisfitSurfaceSummary.Var_Reduc_List;

X_GridForGeogMisfitSurface = [min(Lonstore)-0.5:0.2:max(Lonstore)+0.5];
Y_GridForGeogMisfitSurface = [min(Latstore)-0.5:0.2:max(Latstore)+0.5];
[XXGRD,YYGRD] = meshgrid(X_GridForGeogMisfitSurface,Y_GridForGeogMisfitSurface);

junk2 = figure(149)
nrows = ceil(sqrt(length(EVIDLIST)));
ncols =nrows;
ax = subplot_custom_make(junk2,nrows,ncols,[0.05],[0.05],[0.05 0.95],[0.05 0.95]);

for eventnum = 1:length(RMSList)
    current_RMS = RMSList(eventnum);
    disp(['Completed ' num2str(100*eventnum/length(RMSList)) '%'])
    EVID = EVIDLIST{eventnum}; current_misfit = StoreAllMisfits_L1(:,eventnum);
    FigureName = [MisfitSurfaceFigs_CurrentPeriod 'EVID_' EVID '_MisfitXSections.png'];

    [minmis,mindx] = min(current_misfit);
    % get the model params corresponding to this misfit 
    BestLon = Lonstore(mindx); BestLat = Latstore(mindx);
    BestWidth = Widthstore(mindx); BestTau = Taustore(mindx);
    BestVR = VR_Store(mindx);
     titlestr = {['Event ID ' EVID, ', Best-Fit \tau = ' num2str(BestTau) 's, Best-Fit Width = ' ...
         num2str(BestWidth) 'km'],['Best-Fit Long/Lat; ' num2str(BestLon) '/' num2str(BestLat) '°, Mean(|Observed-Predicted|) = ' num2str(minmis) '°']};

    % Get the 'surfaces' corresponding to this misfit
    VaryingLonLat_Indices = find(Widthstore == BestWidth & ...
        Taustore == BestTau);
    VaryingLons_BestSection = Lonstore(VaryingLonLat_Indices);
    VaryingLats_BestSection = Latstore(VaryingLonLat_Indices);
    L1Misfit_VaryingLonLats= current_misfit(VaryingLonLat_Indices);


    VaryingWidth_Indices = find(Lonstore == BestLon & ...
        Taustore == BestTau & Latstore == BestLat);
    VaryingWidths_BestSection = Widthstore(VaryingWidth_Indices);
    L1Misfit_VaryingWidths = current_misfit(VaryingWidth_Indices);


    VaryingTau_Indices =  find(Lonstore == BestLon & ...
        Widthstore == BestWidth & Latstore == BestLat);
    VaryingTaus_BestSection = Taustore(VaryingTau_Indices);
    L1Misfit_VaryingTaus = current_misfit(VaryingTau_Indices);

    % Now, load the files corresponding to this event 
    % and get the model corresponding to the lowest misfit
    ModelSpaceInfo  = [ ModelStorageFolder...
    'ModelSpaceSearch_Store_EVID' EVID '.mat'];
    load(ModelSpaceInfo)    
    BestFittingModel = Output_Store.ModelSpaceSearch_Store(mindx,:);
    xgrid = Output_Store.xgrid; 
    ygrid = Output_Store.ygrid;
    % Load the observations.    
    RawObs_Fname = [ObservationsDir_CurrentPeriod EVID num2str(Period) 's_elon_elat_lon_lat_phidev_phigc'];
    AngObsInfo =  load(RawObs_Fname,'-ascii');
    Elon = AngObsInfo(1,1); Elat = AngObsInfo(1,2);
    Obs_Xgrid = AngObsInfo(:,3); Obs_Ygrid = AngObsInfo(:,4); AngleResid = AngObsInfo(:,5);
     RawStationData_Fname = [ObservationsDir_CurrentPeriod EVID num2str(Period) 's_slon_slat_stt'];
    RawStnInfo =  load(RawStationData_Fname,'-ascii');  
    RawSlon = RawStnInfo(:,1);     RawSlat = RawStnInfo(:,2);     RawStt = RawStnInfo(:,3);



L1Gridded = griddata(VaryingLons_BestSection,VaryingLats_BestSection,L1Misfit_VaryingLonLats,XXGRD,YYGRD);
contourf(ax(eventnum),XXGRD,YYGRD,L1Gridded,'edgecolor','none')

hold(ax(eventnum),"on")
plot(ax(eventnum),coastlon,coastlat,'color','k','LineWidth',2)
text(ax(eventnum),-162, 18,EVID,'fontsize',16)
   
% Quote and display QC Parameters
text(ax(eventnum),-162, 23,['RMS = ' num2str(round(current_RMS,2)) 's'],'fontsize',16)
  text(ax(eventnum),-162, 20,['N = ' num2str(round(length(Obs_Xgrid),2)) ],'fontsize',16)

%box on
set(ax(eventnum),'linewidth',2,'box','on','XTickLabel',[],'YTickLabel',[])
   ax(eventnum).XLim=[-163 -154];
     ax(eventnum).YLim=[17 25];
ax(eventnum).CLim = [0 7.5];
colormap(ax(eventnum),turbo)
end

%%% 


set(gcf,'position',[1637 -77 1134 956])
currnum = eventnum;

for eventnum = currnum+1:1:nrows*ncols

ax(eventnum).Color = 'none'
ax(eventnum).Box = 'off'
ax(eventnum).XTick = []
ax(eventnum).YTick = []
ax(eventnum).XColor = 'none'
ax(eventnum).YColor = 'none'
ax(eventnum).FontSize = 18

if eventnum == currnum+1
% add colorbar
barbar = colorbar(ax(eventnum),'location','west')
colormap(ax(eventnum),"turbo")
ax(eventnum).CLim = [0 10];
ylabel(barbar,{'Avg. L1 Misfit','(degrees)'}, 'fontsize',18)
end
end
sgtitle([num2str(Period) 's'],'fontsize',20,'fontweight','bold')
saveas(junk2,[FigFolder num2str(Period) '_sAllMisfitSurfaces.png'])
close all
end

