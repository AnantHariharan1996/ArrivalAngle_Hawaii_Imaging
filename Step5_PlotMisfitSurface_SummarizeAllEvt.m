 %% Step 4: Plot Misfit Surfaces
% What does this script do?
% For every period, we stack the misfit for every event.
% Make plots showing the stacked misfit surfaces
% Plot ALL the best-fit plume Models (inc. width) on a map. 
% Color them by variance reduction. 


clear; clc; close all; load coastlines; 
HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';
PredictionsDir =  [HomeDir 'Stored_ModelSpacePredictions/'];
ObservationsDir = [HomeDir 'Raw_ArrivalAngleDeviations/'];
SummaryMisfitDir = [PredictionsDir 'SummaryMisfitStore/'];
FigFolder =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/Stored_ModelSpacePredictions/SummaryFigures/';
Periodlist = 66.6667;
Pcounter = 0;
for  Period=Periodlist
    Pcounter = Pcounter +1;
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

X_GridForGeogMisfitSurface = [min(Lonstore)-0.5:0.2:max(Lonstore)+0.5];
Y_GridForGeogMisfitSurface = [min(Latstore)-0.5:0.2:max(Latstore)+0.5];
[XXGRD,YYGRD] = meshgrid(X_GridForGeogMisfitSurface,Y_GridForGeogMisfitSurface);

for eventnum = 1:length(RMSList)

    disp(['Completed ' num2str(100*eventnum/length(RMSList)) '%'])
    EVID = EVIDLIST{eventnum}; current_misfit = StoreAllMisfits_L1(:,eventnum);
    FigureName = [MisfitSurfaceFigs_CurrentPeriod 'EVID_' EVID '_MisfitXSections.png'];

    [minmis,mindx] = min(current_misfit);
    % get the model params corresponding to this misfit 
    BestLon = Lonstore(mindx); BestLat = Latstore(mindx);
    BestWidth = Widthstore(mindx); BestTau = Taustore(mindx);
     titlestr = {['Event ID ' EVID, ', Best-Fit \tau = ' num2str(BestTau) 's, Best-Fit Width = ' ...
         num2str(BestWidth) 'km'],['Best-Fit Long/Lat; ' num2str(BestLon) '/' num2str(BestLat) '°' ]};

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
    
    % Map showing the best-fit plume locations
    figure(5000+Pcounter)
       scatter(BestLon, BestLat,250,BestTau,'filled')
       hold on


    figure(8000+Pcounter)
       scatter(BestLon, BestLat,250,minmis,'filled')
hold on

end
zzz(:,1) = MisfitSurfaceSummary.Lonstore;
zzz(:,2) = MisfitSurfaceSummary.Latstore;

    figure(5000+Pcounter)

box on
    plot(coastlon,coastlat,'linewidth',2,'color','k')
    hold on
barbar=colorbar;
ylabel(barbar,'timedelay (s)')

[uniquerows,idx] = unique(zzz,'rows');
       scatter(uniquerows(:,1),uniquerows(:,2),1000,[0 0 1],'s','filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.05)
    set(gca,'linewidth',2,'fontsize',17,'layer','top','fontweight','bold')
       xlim([-163 -153])
       ylim([16 26])

       hold on
colormap(turbo)
saveas(figure(5000+Pcounter),[FigFolder num2str(Period) 's_alllocssummary_WithTau.png'])



    figure(8000+Pcounter)

box on
    plot(coastlon,coastlat,'linewidth',2,'color','k')
    hold on
barbar=colorbar;
ylabel(barbar,'Misfit (deg)')

[uniquerows,idx] = unique(zzz,'rows');
       scatter(uniquerows(:,1),uniquerows(:,2),1000,[0 0 1],'s','filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.05)
    set(gca,'linewidth',2,'fontsize',17,'layer','top','fontweight','bold')
       xlim([-163 -153])
       ylim([16 26])

       hold on
colormap(turbo)
saveas(figure(8000+Pcounter),[FigFolder num2str(Period) 's_alllocssummary_withMisfit.png'])


% Stack the Misfit. 
L1_MisfitSurfaceSummary_stacked =sum(StoreAllMisfits_L1');
L2_MisfitSurfaceSummary_stacked =sum(StoreAllMisfits_L2');

% get cross-sections for plotting
[minmis,mindx_L1] = min(L1_MisfitSurfaceSummary_stacked);
BestLon_L1 = Lonstore(mindx_L1); BestLat_L1 = Latstore(mindx_L1);
BestWidth_L1 = Widthstore(mindx_L1); BestTau_L1 = Taustore(mindx_L1);

[minmis,mindx_L2] = min(L2_MisfitSurfaceSummary_stacked);
BestLon_L2 = Lonstore(mindx_L2); BestLat_L2 = Latstore(mindx_L2);
BestWidth_L2 = Widthstore(mindx_L2); BestTau_L2 = Taustore(mindx_L2);

% 
L1_VaryingLonLat_Indices = find(Widthstore == BestWidth_L1 & ...
    Taustore == BestTau_L1);
L1_VaryingLons_BestSection = Lonstore(L1_VaryingLonLat_Indices);
L1_VaryingLats_BestSection = Latstore(L1_VaryingLonLat_Indices);
L1Misfit_VaryingLonLats = L1_MisfitSurfaceSummary_stacked(L1_VaryingLonLat_Indices);
% 
L2_VaryingLonLat_Indices = find(Widthstore == BestWidth_L2 & ...
    Taustore == BestTau_L2);
L2_VaryingLons_BestSection = Lonstore(L2_VaryingLonLat_Indices);
L2_VaryingLats_BestSection = Latstore(L2_VaryingLonLat_Indices);
L2Misfit_VaryingLonLats = L2_MisfitSurfaceSummary_stacked(L2_VaryingLonLat_Indices);

% 
L1Misfit_VaryingTau_Indices = find(Widthstore == BestWidth_L1 & ...
    Lonstore == BestLon_L1 & Latstore == BestLat_L1);
L1Misfit_VaryingWidth_Indices = find(Taustore == BestTau_L1 & ...
    Lonstore == BestLon_L1 & Latstore == BestLat_L1);

L2Misfit_VaryingTau_Indices = find(Widthstore == BestWidth_L2 & ...
    Lonstore == BestLon_L2 & Latstore == BestLat_L2);
L2Misfit_VaryingWidth_Indices = find(Taustore == BestTau_L2 & ...
    Lonstore == BestLon_L2 & Latstore == BestLat_L2);

L1_VaryingTau_BestSection = Taustore(L1Misfit_VaryingTau_Indices);
L1Misfit_VaryingTau = L1_MisfitSurfaceSummary_stacked(L1Misfit_VaryingTau_Indices);
L1_VaryingWidth_BestSection = Widthstore(L1Misfit_VaryingWidth_Indices);
L1Misfit_VaryingWidth = L1_MisfitSurfaceSummary_stacked(L1Misfit_VaryingWidth_Indices);


L2_VaryingTau_BestSection = Taustore(L2Misfit_VaryingTau_Indices);
L2Misfit_VaryingTau = L2_MisfitSurfaceSummary_stacked(L2Misfit_VaryingTau_Indices);
L2_VaryingWidth_BestSection = Widthstore(L2Misfit_VaryingWidth_Indices);
L2Misfit_VaryingWidth = L2_MisfitSurfaceSummary_stacked(L2Misfit_VaryingWidth_Indices);



%% 

L1Gridded = griddata(L1_VaryingLons_BestSection,L1_VaryingLats_BestSection,L1Misfit_VaryingLonLats./max(L1Misfit_VaryingLonLats),XXGRD,YYGRD);
L2Gridded = griddata(L2_VaryingLons_BestSection,L2_VaryingLats_BestSection,L2Misfit_VaryingLonLats./max(L2Misfit_VaryingLonLats),XXGRD,YYGRD);

figure(9000+Pcounter)
subplot(2,2,1)
    plot(coastlon,coastlat,'linewidth',2,'color','k')
    hold on
contourf(XXGRD,YYGRD,L1Gridded,50,'edgecolor','none')
barbar=colorbar;
ylabel(barbar,'Normalized L1 Misfit')
box on
set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
title({'Stacked L1 Misfit Surface: Varying Position',['\tau = ' num2str(BestTau_L1), 's, Width = ' num2str(BestWidth_L1) ' km']})
       xlim([-163 -154])
       ylim([17 25])
plot(coastlon,coastlat,'linewidth',3,'color','k')
Cmap2Use= '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/UsefulFunctions/roma.cpt';  

 cptcmap(Cmap2Use,'ncol',20);


subplot(2,2,2)
    plot(coastlon,coastlat,'linewidth',2,'color','k')
hold on
title({'Stacked L2 Misfit Surface: Varying Position',['\tau = ' num2str(BestTau_L2), 's, Width = ' num2str(BestWidth_L2) ' km']})
contourf(XXGRD,YYGRD,L2Gridded,50,'edgecolor','none')
barbar=colorbar;
ylabel(barbar,'Normalized L2 Misfit')
set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
box on
plot(coastlon,coastlat,'linewidth',3,'color','k')
       xlim([-163 -154])
       ylim([17 25])
Cmap2Use= '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/UsefulFunctions/roma.cpt';  

 cptcmap(Cmap2Use,'ncol',20);

       set(gcf,'position', [37 43 1041 791])


%% Now plot the misfit surfaces for Tau and for Width
subplot(2,2,3)
yyaxis left 
    plot(L1_VaryingTau_BestSection,L1Misfit_VaryingTau,'linewidth',2)
    ylabel('L1 Misfit')
yyaxis right
    plot(L2_VaryingTau_BestSection,L2Misfit_VaryingTau,'linewidth',2)
xlabel('\tau (s)')
    ylabel('L2 Misfit')
    box on
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)

    title('Stacked Misfit Surface: Varying Time Delay')


subplot(2,2,4)

yyaxis left 
    plot(L1_VaryingWidth_BestSection,L1Misfit_VaryingWidth,'linewidth',2)
    ylabel('L1 Misfit')
yyaxis right
    plot(L2_VaryingWidth_BestSection,L2Misfit_VaryingWidth,'linewidth',2)
xlabel('Width (km)')
    ylabel('L2 Misfit')
    box on
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
    title('Stacked Misfit Surfaces: Varying Width')
saveas(figure(9000+Pcounter),[FigFolder  num2str(Period) 's_StackedMisfitSurface.png'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

