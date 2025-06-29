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
Periodlist = 50;
Pcounter = 0;
PrctileThresh = 2;
addpath(genpath(pwd))
TauBinsforHist = [0:5:40];
WidthBinsforHist = [0:50:400];
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
StoreAllMisfits_L1_Weighted = MisfitSurfaceSummary.StoreAllMisfits_L1_Weighted;
StoreAllMisfits_L2_Weighted = MisfitSurfaceSummary.StoreAllMisfits_L2_Weighted;
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


zzz(:,1) = MisfitSurfaceSummary.Lonstore;
zzz(:,2) = MisfitSurfaceSummary.Latstore;


[uniquerows,idx] = unique(zzz,'rows');



% Stack the Misfit. 
L1_MisfitSurfaceSummary_stacked =mean(StoreAllMisfits_L1');
L2_MisfitSurfaceSummary_stacked =mean(StoreAllMisfits_L2');

L1_MisfitSurfaceSummary_stacked_weighted = mean(StoreAllMisfits_L1_Weighted');
L2_MisfitSurfaceSummary_stacked_weighted= mean(StoreAllMisfits_L2_Weighted');


% get cross-sections for plotting
[minmisL1,mindx_L1] = min(L1_MisfitSurfaceSummary_stacked);
BestLon_L1 = Lonstore(mindx_L1); BestLat_L1 = Latstore(mindx_L1);
BestWidth_L1 = Widthstore(mindx_L1); BestTau_L1 = Taustore(mindx_L1);

[minmisL2,mindx_L2] = min(L2_MisfitSurfaceSummary_stacked);
BestLon_L2 = Lonstore(mindx_L2); BestLat_L2 = Latstore(mindx_L2);
BestWidth_L2 = Widthstore(mindx_L2); BestTau_L2 = Taustore(mindx_L2);


[minmisL1_W,mindx_L1_W] = min(L1_MisfitSurfaceSummary_stacked_weighted);
BestLon_L1_W = Lonstore(mindx_L1_W); 
BestLat_L1_W = Latstore(mindx_L1_W);
BestWidth_L1_W = Widthstore(mindx_L1_W); 
BestTau_L1_W = Taustore(mindx_L1_W);


[minmisL2_W,mindx_L2_W] = min(L2_MisfitSurfaceSummary_stacked_weighted);
BestLon_L2_W = Lonstore(mindx_L2_W); 
BestLat_L2_W = Latstore(mindx_L2_W);
BestWidth_L2_W = Widthstore(mindx_L2_W); 
BestTau_L2_W = Taustore(mindx_L2_W);


% get transects along lon and lat
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

% now get the transects along lon and lat for the weighted models
L1_VaryingLonLat_Indices_W = find(Widthstore == BestWidth_L1_W & ...
    Taustore == BestTau_L1_W);
L1_VaryingLons_BestSection_W = Lonstore(L1_VaryingLonLat_Indices_W);
L1_VaryingLats_BestSection_W = Latstore(L1_VaryingLonLat_Indices_W);
L1Misfit_VaryingLonLats_W = L1_MisfitSurfaceSummary_stacked(L1_VaryingLonLat_Indices_W);
% 

L2_VaryingLonLat_Indices_W = find(Widthstore == BestWidth_L2_W & ...
    Taustore == BestTau_L2_W);
L2_VaryingLons_BestSection_W = Lonstore(L2_VaryingLonLat_Indices_W);
L2_VaryingLats_BestSection_W = Latstore(L2_VaryingLonLat_Indices_W);
L2Misfit_VaryingLonLats_W = L2_MisfitSurfaceSummary_stacked(L2_VaryingLonLat_Indices_W);
% 



L1Misfit_VaryingTau_Indices = find(Widthstore == BestWidth_L1 & ...
    Lonstore == BestLon_L1 & Latstore == BestLat_L1);
L1Misfit_VaryingWidth_Indices = find(Taustore == BestTau_L1 & ...
    Lonstore == BestLon_L1 & Latstore == BestLat_L1);

L1Misfit_VaryingTau_Indices_W = find(Widthstore == BestWidth_L1_W & ...
    Lonstore == BestLon_L1_W & Latstore == BestLat_L1_W);
L1Misfit_VaryingWidth_Indices_W = find(Taustore == BestTau_L1_W & ...
    Lonstore == BestLon_L1_W & Latstore == BestLat_L1_W);

L2Misfit_VaryingTau_Indices = find(Widthstore == BestWidth_L2 & ...
    Lonstore == BestLon_L2 & Latstore == BestLat_L2);
L2Misfit_VaryingWidth_Indices = find(Taustore == BestTau_L2 & ...
    Lonstore == BestLon_L2 & Latstore == BestLat_L2);

L2Misfit_VaryingTau_Indices_W = find(Widthstore == BestWidth_L2_W & ...
    Lonstore == BestLon_L2_W & Latstore == BestLat_L2_W);
L2Misfit_VaryingWidth_Indices_W = find(Taustore == BestTau_L2_W & ...
    Lonstore == BestLon_L2_W & Latstore == BestLat_L2_W);


L1_VaryingTau_BestSection = Taustore(L1Misfit_VaryingTau_Indices);
L1Misfit_VaryingTau = L1_MisfitSurfaceSummary_stacked(L1Misfit_VaryingTau_Indices);
L1_VaryingWidth_BestSection = Widthstore(L1Misfit_VaryingWidth_Indices);
L1Misfit_VaryingWidth = L1_MisfitSurfaceSummary_stacked(L1Misfit_VaryingWidth_Indices);

L2_VaryingTau_BestSection = Taustore(L2Misfit_VaryingTau_Indices);
L2Misfit_VaryingTau = L2_MisfitSurfaceSummary_stacked(L2Misfit_VaryingTau_Indices);
L2_VaryingWidth_BestSection = Widthstore(L2Misfit_VaryingWidth_Indices);
L2Misfit_VaryingWidth = L2_MisfitSurfaceSummary_stacked(L2Misfit_VaryingWidth_Indices);




%% Extract distributions for plotting.
MinMisfit_L1 = min(L1_MisfitSurfaceSummary_stacked);
MinMisfit_L2 = min(L2_MisfitSurfaceSummary_stacked);
MinMisfit_L1_W = min(L1_MisfitSurfaceSummary_stacked_weighted);
MinMisfit_L2_W= min(L2_MisfitSurfaceSummary_stacked_weighted);

Low_L1_MisfitThreshold = prctile(L1_MisfitSurfaceSummary_stacked,PrctileThresh);
Low_L2_MisfitThreshold = prctile(L2_MisfitSurfaceSummary_stacked,PrctileThresh);
Low_L1_MisfitThreshold_W = prctile(L1_MisfitSurfaceSummary_stacked_weighted,PrctileThresh);
Low_L2_MisfitThreshold_W = prctile(L2_MisfitSurfaceSummary_stacked_weighted,PrctileThresh);

BestFit_idx_L1 = find(L1_MisfitSurfaceSummary_stacked < Low_L1_MisfitThreshold);
BestFit_idx_L2 = find(L2_MisfitSurfaceSummary_stacked < Low_L2_MisfitThreshold);
BestFit_idx_L1_W = find(L1_MisfitSurfaceSummary_stacked_weighted < Low_L1_MisfitThreshold_W);
BestFit_idx_L2_W = find(L2_MisfitSurfaceSummary_stacked_weighted < Low_L2_MisfitThreshold_W);

%% Tau distros
TopPercentile_Tau_L1 = Taustore(BestFit_idx_L1);
TopPercentile_Tau_L2 = Taustore(BestFit_idx_L2);
TopPercentile_Tau_L1_W = Taustore(BestFit_idx_L1_W);
TopPercentile_Tau_L2_W = Taustore(BestFit_idx_L2_W);

%% Width distros
TopPercentile_Width_L1 = Widthstore(BestFit_idx_L1);
TopPercentile_Width_L2 = Widthstore(BestFit_idx_L2);
TopPercentile_Width_L1_W = Widthstore(BestFit_idx_L1_W);
TopPercentile_Width_L2_W = Widthstore(BestFit_idx_L2_W);

%% Position distros
TopPercentile_Lon_L1 = Lonstore(BestFit_idx_L1);
TopPercentile_Lon_L2 = Lonstore(BestFit_idx_L2);
TopPercentile_Lon_L1_W = Lonstore(BestFit_idx_L1_W);
TopPercentile_Lon_L2_W = Lonstore(BestFit_idx_L2_W);

TopPercentile_Lat_L1 = Latstore(BestFit_idx_L1);
TopPercentile_Lat_L2 = Latstore(BestFit_idx_L2);
TopPercentile_Lat_L1_W = Latstore(BestFit_idx_L1_W);
TopPercentile_Lat_L2_W = Latstore(BestFit_idx_L2_W);

L1GoodLocs(:,1) = TopPercentile_Lon_L1;
L1GoodLocs(:,2) = TopPercentile_Lat_L1;

L2GoodLocs(:,1) = TopPercentile_Lon_L2;
L2GoodLocs(:,2) = TopPercentile_Lat_L2;

L1GoodLocs_W(:,1) = TopPercentile_Lon_L1_W;
L1GoodLocs_W(:,2) = TopPercentile_Lat_L1_W;

L2GoodLocs_W(:,1) = TopPercentile_Lon_L2_W;
L2GoodLocs_W(:,2) = TopPercentile_Lat_L2_W;

[unique_L1,unique_L1dx] = unique(L1GoodLocs,'rows');
[unique_L2,unique_L2dx] = unique(L2GoodLocs,'rows');
[unique_L1_W,unique_L1dx_W] = unique(L1GoodLocs_W,'rows');
[unique_L2_W,unique_L2dx_W] = unique(L2GoodLocs_W,'rows');

% Loop over every unique point and get the number of matches. 
for uniquelocnum = 1:length(unique_L1(:,1))
    current_lon = unique_L1(uniquelocnum,1);
    current_lat = unique_L1(uniquelocnum,2);
    numvals_L1(uniquelocnum) = length(find( ...
        TopPercentile_Lon_L1 == current_lon & TopPercentile_Lat_L1 == current_lat));
end

% Loop over every unique point and get the number of matches. 
for uniquelocnum = 1:length(unique_L2(:,1))
    current_lon = unique_L2(uniquelocnum,1);
    current_lat = unique_L2(uniquelocnum,2);
    numvals_L2(uniquelocnum) = length(find( ...
        TopPercentile_Lon_L2 == current_lon & TopPercentile_Lat_L2 == current_lat));
end

% Loop over every unique point and get the number of matches. 
for uniquelocnum = 1:length(unique_L1_W(:,1))

    current_lon = unique_L1_W(uniquelocnum,1);
    current_lat = unique_L1_W(uniquelocnum,2);
    numvals_L1_W(uniquelocnum) = length(find( ...
    TopPercentile_Lon_L1_W == current_lon & TopPercentile_Lat_L2_W == current_lat));

end


% Loop over every unique point and get the number of matches. 
for uniquelocnum = 1:length(unique_L2_W(:,1))

    current_lon = unique_L2_W(uniquelocnum,1);
    current_lat = unique_L2_W(uniquelocnum,2);
    numvals_L2_W(uniquelocnum) = length(find( ...
    TopPercentile_Lon_L2_W == current_lon & TopPercentile_Lat_L2_W == current_lat));

end

%% 

L1Gridded = griddata(L1_VaryingLons_BestSection,L1_VaryingLats_BestSection,L1Misfit_VaryingLonLats,XXGRD,YYGRD);
L2Gridded = griddata(L2_VaryingLons_BestSection,L2_VaryingLats_BestSection,L2Misfit_VaryingLonLats,XXGRD,YYGRD);




%%% ALL PLOTTING BELOW HERE!
figure(9000+Pcounter)
subplot(2,3,1)

    plot(coastlon,coastlat,'linewidth',2,'color','k')
hold on
title({'Stacked L2 Misfit Surface: Varying Position',['\tau = ' num2str(BestTau_L2), 's, Width = ' num2str(BestWidth_L2) ' km']})
contourf(XXGRD,YYGRD,L2Gridded,50,'edgecolor','none')
barbar=colorbar;
ylabel(barbar,'L2 Misfit')
set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
box on
plot(coastlon,coastlat,'linewidth',3,'color','k')
       xlim([-163 -154])
       ylim([17 25])
Cmap2Use= '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/UsefulFunctions/roma.cpt';  

 cptcmap(Cmap2Use,'ncol',20);

       set(gcf,'position', [37 43 1041 791])


%% Now plot the misfit surfaces for Tau and for Width
subplot(2,3,2)
yyaxis left 
    plot(L1_VaryingTau_BestSection,L1Misfit_VaryingTau,'-o','linewidth',2)
    ylabel('L1 Misfit')
yyaxis right
    plot(L2_VaryingTau_BestSection,L2Misfit_VaryingTau,'-o','linewidth',2)
xlabel('\tau (s)')
    ylabel('L2 Misfit')
    box on
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
grid on;
    title('Stacked Misfit Surface: Varying Time Delay')


subplot(2,3,3)

yyaxis left 
    plot(L1_VaryingWidth_BestSection,L1Misfit_VaryingWidth,'-o','linewidth',2)
    ylabel('L1 Misfit')
yyaxis right
    plot(L2_VaryingWidth_BestSection,L2Misfit_VaryingWidth,'-o','linewidth',2)
xlabel('Width (km)')
    ylabel('L2 Misfit')
    box on
    grid on;
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
    title('Stacked Misfit Surfaces: Varying Width')
    set(figure(9000+Pcounter),'position', [1600 54 1576 789])


subplot(2,3,5)
histogram(TopPercentile_Tau_L1,TauBinsforHist)
hold on
histogram(TopPercentile_Tau_L2,TauBinsforHist)
% histogram(TopPercentile_Tau_L1_W,TauBinsforHist)
% histogram(TopPercentile_Tau_L1_W,TauBinsforHist)
title(['Top ' num2str(PrctileThresh) ' Percentile of Models'])
xlabel('\tau (s)')
box on
grid on;
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
    legend('L1','L2')

subplot(2,3,6)
histogram(TopPercentile_Width_L1,WidthBinsforHist)
hold on
histogram(TopPercentile_Width_L2,WidthBinsforHist)
% histogram(TopPercentile_Width_L1_W,WidthBinsforHist)
% histogram(TopPercentile_Width_L1_W,WidthBinsforHist)
title(['Top ' num2str(PrctileThresh) ' Percentile of Models'])
xlabel('Width (km)')
box on
grid on;
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
    legend('L1','L2','location','northwest')

    %%% Make the 2D histogram
 
subplot(2,3,4)
scatter(unique_L2(:,1),unique_L2(:,2),150,numvals_L2,'filled')
hold on
plot(coastlon,coastlat,'linewidth',3,'color','k')
       xlim([-163 -154])
       ylim([17 25])
barbar=colorbar;
box on
ylabel(barbar,'N')
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
title(['Top ' num2str(PrctileThresh) ' Percentile of Models: L2'])
saveas(figure(9000+Pcounter),[FigFolder  num2str(Period) 's_StackedMisfitSurface.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Now, repeat the above plot, but this time for the weighted misfit%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(900+Pcounter)

subplot(2,3,1)

subplot(2,3,2)

subplot(2,3,3)

subplot(2,3,4)

subplot(2,3,5)

histogram(TopPercentile_Tau_L1_W,TauBinsforHist)
hold om
histogram(TopPercentile_Tau_L2_W,TauBinsforHist)
title(['Top ' num2str(PrctileThresh) ' Percentile of Models'])
xlabel('\tau (s)')
box on
grid on;
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
    legend('L1','L2')

subplot(2,3,6)

histogram(TopPercentile_Width_L1_W,WidthBinsforHist)
hold on
histogram(TopPercentile_Width_L2_W,WidthBinsforHist)
title(['Top ' num2str(PrctileThresh) ' Percentile of Models'])
xlabel('Width (km)')
box on
grid on;
set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
legend('L1','L2','location','northwest')

saveas(figure(9000+Pcounter),[FigFolder  num2str(Period) 's_StackedMisfitSurfaceWeighted.png'])

end

