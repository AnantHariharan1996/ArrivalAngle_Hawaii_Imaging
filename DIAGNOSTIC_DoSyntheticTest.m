%% Do synthetic test
% For an assumed synthetic diffractor location, time delay, etc. 
% See how well it could be recovered, given the geometry we have for our
% problem. 


clear; clc; close all;
SyntheticTestFolder =     '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/Stored_ModelSpacePredictions/SyntheticTests/';
HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';
Storagefolder = [HomeDir 'Stored_ModelSpacePredictions/']
% Step 0: What period?
Period  = 50;
Storagefolder = [Storagefolder num2str(Period) 's/'];
% some settings
[phvel]=prem_dispersion(1/Period);
cglb = phvel;
spacing=0.25;
lambda = Period*cglb;
load coastlines 
NoiseFac = 5;
% 

PlotTrueObservations = 1; 
if PlotTrueObservations
junk2 = figure(90210)
nrows = 8; ncols = 7;
ax = subplot_custom_make(junk2,nrows,ncols,[0.05],[0.05],[0.05 0.95],[0.05 0.95]);
end


NewFolder  = [HomeDir 'Raw_ArrivalAngleDeviations/' num2str(Period) 's/'];

% Step 1: Define 'true' diffractor parameters. Make a plot showing this. 
True_Scatterer_Lon = [-160];
True_Scatterer_Lat= [22];
True_Scatterer_Width = [100];
True_Scatterer_Tau = [30];


Event_Files=  dir([NewFolder '*_elon_elat_lon_lat_phidev_phigc']);
Predictions_List = dir([Storagefolder 'ModelSpaceSearch_Store_EVID*'])


for EventNum = 1:length(Predictions_List)
disp(['Event: ' num2str(EventNum)])
    % Step 2: For every earthquake in the dataset, run the forward modeling
% code and generate 'observations' 
% AT THE LOCATIONS WHERE WE HAVE OUR OBSERVATIONS


       Current_Event_File = Event_Files(EventNum).name;
       % Get Event ID and geometry
       CurrID =extractBefore(Current_Event_File,[num2str(Period) 's']);
       CurrID_num =  str2num(CurrID); CurrInfo = load([NewFolder Current_Event_File],'-ascii');
       EvLon = CurrInfo(1,1); EvLat = CurrInfo(1,2);
       CurrRawDataInfo = load([NewFolder CurrID  num2str(Period) 's' '_slon_slat_stt'],'-ascii');
       FullFname_Current_Event_File = [NewFolder Current_Event_File];
        
       ObsInfo = load(FullFname_Current_Event_File,'-ascii');
       Obs_Lon = ObsInfo(:,3);  Obs_Lat = ObsInfo(:,4); Obs_Delta = ObsInfo(:,5);
        
       % Generate the Grid used for the predictions
         Ref_XGrid = min(Obs_Lon)-0.25:spacing:max(Obs_Lon)+0.25;
         Ref_YGrid = min(Obs_Lat)-0.25:spacing:max(Obs_Lat)+0.25;
          [Ref_XGrid,Ref_YGrid] = meshgrid(Ref_XGrid,Ref_YGrid);

         % Running forward model! %

         scattererlat  = True_Scatterer_Lat;
         scattererlon = True_Scatterer_Lon;
         current_timelag = True_Scatterer_Tau;
         current_width = True_Scatterer_Width;
        [delta,xgrid,ygrid,ttime_field_perturbed,tau,ttime_field_noscatter,angle,angle_nodiff] ...
            = Get_Arrival_Angle_Residual_GaussianBeam(Period, EvLat,...
            EvLon, scattererlat,scattererlon,current_timelag,current_width,...
            Ref_YGrid(:),Ref_XGrid(:),cglb,spacing);


        InterpedAAsOnObs = griddata(xgrid,ygrid,delta,Obs_Lon,Obs_Lat);
        Noise = rand(size(InterpedAAsOnObs))-0.5;
        InterpedAAsOnObs = InterpedAAsOnObs+NoiseFac.*Noise;
        Interped_AAs = griddata(xgrid,ygrid,delta,Ref_XGrid,Ref_YGrid);

        if PlotTrueObservations
            % to verify things are reasonable, let's plot the forward
            % model.    

        % contourf(ax(EventNum),Ref_XGrid,Ref_YGrid,interped_AAs,100,'edgecolor','none')
        scatter(ax(EventNum),Obs_Lon,Obs_Lat,5,InterpedAAsOnObs,'filled')
        hold(ax(EventNum),"on")
        plot(ax(EventNum),coastlon,coastlat,'color','k','LineWidth',2)
        % text(ax(EventNum),-162, 18,CurrID,'fontsize',16)
           scatter(ax(EventNum),True_Scatterer_Lon,True_Scatterer_Lat,250,[1 0 0],'filled')
        ax(EventNum).CLim = [-5 5]; 
        ax(EventNum).XLim=[-163 -150];
        ax(EventNum).YLim=[17 25];
        colormap(ax(EventNum),turbo(50))
        ax(EventNum).XTick = [];
        ax(EventNum).YTick = [];

        end

    % Step 3: Loop over EVERY model parameter for this observation and
    % calculate the misfit


        % Now, load the database of predictions for this event
        Predictions_Name =[Storagefolder 'ModelSpaceSearch_Store_EVID'  CurrID '.mat'];

        load(Predictions_Name)

            Lonstore = Output_Store.Lonstore;
          Latstore = Output_Store.Latstore;
          Widthstore = Output_Store.Widthstore;
          Taustore = Output_Store.Taustore;

X_GridForGeogMisfitSurface = [min(Lonstore)-0.5:0.2:max(Lonstore)+0.5];
Y_GridForGeogMisfitSurface = [min(Latstore)-0.5:0.2:max(Latstore)+0.5];
[XXGRD,YYGRD] = meshgrid(X_GridForGeogMisfitSurface,Y_GridForGeogMisfitSurface);

    L1_MisfitStore = 999999.*ones(size(Taustore)); 
           L2_MisfitStore = L1_MisfitStore;
           Var_Reduc_List = L1_MisfitStore;
           
          for modelspace_num = 1:length(Taustore)
                
                current_predictions = Output_Store.ModelSpaceSearch_Store(modelspace_num,:);
                % Get misfit between 'observations' and predictions
                 Interped_AA_PREDICTIONS = griddata(Output_Store.xgrid,Output_Store.ygrid,current_predictions,...
                      Obs_Lon,Obs_Lat);
                Difference = InterpedAAsOnObs-Interped_AA_PREDICTIONS;
                Mean_L1Misfit = mean(abs(Difference)); % currently unweighted (but still avg) misfit
                Mean_L2Misfit = mean(abs(Difference.^2)); % currently unweighted (but still avg) misfit
                L1_MisfitStore(modelspace_num) = Mean_L1Misfit;
                L2_MisfitStore(modelspace_num) = Mean_L2Misfit;         
          end

          StoreAllMisfits_L1(:,EventNum) = L1_MisfitStore;
          StoreAllMisfits_L2(:,EventNum) = L2_MisfitStore;
          

end

% clean up the figure that shows all the synthetic observations




% Misfits are now stored in massive matrices. 

% Stack the Misfit. 
L1_MisfitSurfaceSummary_stacked =sum(StoreAllMisfits_L1');
L2_MisfitSurfaceSummary_stacked =sum(StoreAllMisfits_L2');

%% Recover the model parameters that correspond to the best-fitting models 

% Get the top 5% best-fitting models
Low_L1_MisfitThreshold = prctile(L1_MisfitSurfaceSummary_stacked,1);
BestFittingModelDX = ...
    find(L1_MisfitSurfaceSummary_stacked < Low_L1_MisfitThreshold);

Best_Fitting_Taus = Taustore(BestFittingModelDX);
Best_Fitting_Widths = Widthstore(BestFittingModelDX)
Best_Fitting_Lon_L1 = Lonstore(BestFittingModelDX); 
Best_Fitting_Lat_L1 = Latstore(BestFittingModelDX);

% Generate 2D histogram in position space
TempPositionStore(:,1) = Best_Fitting_Lon_L1;
TempPositionStore(:,2) = Best_Fitting_Lat_L1;

[UniquePositionVec,uniquedx] = unique(TempPositionStore,'rows');
for histogram_val = 1:length(UniquePositionVec)
currpos = UniquePositionVec(histogram_val,:);
currposLon = currpos(1); currposLat = currpos(2);
Numdx = length(find(Best_Fitting_Lon_L1 == currposLon & Best_Fitting_Lat_L1 == currposLat));
BestLonList(histogram_val) = currposLon; 
BestLatList(histogram_val) = currposLat;
NumVals(histogram_val) = (Numdx);

end

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

L1Gridded = griddata(L1_VaryingLons_BestSection,L1_VaryingLats_BestSection,L1Misfit_VaryingLonLats,XXGRD,YYGRD);
L2Gridded = griddata(L2_VaryingLons_BestSection,L2_VaryingLats_BestSection,L2Misfit_VaryingLonLats,XXGRD,YYGRD);


%%%%%%%%%PLOTTING BELOW HERE%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9000+Period)
subplot(2,3,1)
plot(coastlon,coastlat,'linewidth',3,'color','k')
    hold on
[m,p1]=contourf(XXGRD,YYGRD,L1Gridded,50,'edgecolor','none');
barbar=colorbar;
ylabel(barbar,'Normalized L1 Misfit')
box on
plot(coastlon,coastlat,'linewidth',3,'color','k')

set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
title({'Stacked L1 Misfit Surface: Varying Position',['\tau = ' num2str(BestTau_L1), 's, Width = ' num2str(BestWidth_L1) ' km']})
       xlim([-165 -154])
       ylim([17 26])
Cmap2Use= '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/UsefulFunctions/roma.cpt';  

 cptcmap(Cmap2Use,'ncol',20);
p2= scatter(True_Scatterer_Lon,True_Scatterer_Lat,500,[1 0 1],'x','linewidth',5);
legend([p1 p2],'Misfit Slice','True Location')
% True_Scatterer_Lon = [-160];
% True_Scatterer_Lat= [20];


%% Now plot the misfit surfaces for Tau and for Width
subplot(2,3,2)

yyaxis left 
    p3=plot(L1_VaryingTau_BestSection,L1Misfit_VaryingTau,'linewidth',2)
    ylabel('L1 Misfit')
yyaxis right
    p4 = plot(L2_VaryingTau_BestSection,L2Misfit_VaryingTau,'linewidth',2)
xlabel('\tau (s)')
    ylabel('L2 Misfit')
    box on
    hold on
set(gca,'fontsize',18,'fontweight','bold','linewidth',2)

    title('Stacked Misfit Surface: Varying Time Delay')

p5 = plot([True_Scatterer_Tau True_Scatterer_Tau],[0 max(L2Misfit_VaryingTau)],'linewidth',3,'LineStyle','--')
legend([p3 p4 p5],'L1 Misfit Slice','L2 Misfit Slice','True \tau')
subplot(2,3,3)

yyaxis left 
   p6= plot(L1_VaryingWidth_BestSection,L1Misfit_VaryingWidth,'linewidth',2);
    ylabel('L1 Misfit')
yyaxis right
hold on
    p7 = plot(L2_VaryingWidth_BestSection,L2Misfit_VaryingWidth,'linewidth',2);
xlabel('Width (km)')
    ylabel('L2 Misfit')
    box on
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
    title('Stacked Misfit Surfaces: Varying Width')
p8 = plot([True_Scatterer_Width True_Scatterer_Width],[0 max(L2Misfit_VaryingWidth)],'linewidth',3,'LineStyle','--');
legend([p6 p7 p8],'L1 Misfit Slice','L2 Misfit Slice','True \tau')


% L2Misfit_VaryingWidth

subplot(2,3,5)
histogram(Best_Fitting_Taus,50)
hold on
plot([True_Scatterer_Tau True_Scatterer_Tau],[0 25],'linewidth',3,'LineStyle','--')
% legend('95th Percentile Best Model Ensemble','True Model')
title({'Histogram of \tau for','top 1 Percentile Best-Fit Values'})
xlabel('\tau (s)')
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
xlim([min(Taustore) max(Taustore)])

subplot(2,3,6)


histogram(Best_Fitting_Widths,50)
hold on
plot([True_Scatterer_Width True_Scatterer_Width],[0 25],'linewidth',3,'LineStyle','--')
% legend('95th Percentile Best Model Ensemble','True Model')
title({'Histogram of Width for','top 1 Percentile Best-Fit Values'})
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)

subplot(2,3,4)
plot(coastlon,coastlat,'linewidth',3,'color','k')
hold on
scatter(BestLonList,BestLatList,150,NumVals,'filled')
set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
title({'2D Histogram of Locations for','top 1 Percentile Best-Fit Values'})
       xlim([-165 -154])
       ylim([17 26])
    set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
tinybar=colorbar;
ylabel(tinybar,'N(models)')


       set(gcf,'position',  [1600 54 1576 789])

sgtitle(['Synthetic Test Results: Diffractor Location' num2str(True_Scatterer_Lon) '/' num2str(True_Scatterer_Lat) ', \tau = ' num2str(True_Scatterer_Tau) 's, Width = ' num2str(True_Scatterer_Width) ' km'],'fontsize',24,'fontweight','bold')
saveas(figure(9000+Period),[SyntheticTestFolder  num2str(Period) 's_StackedMisfitSurfaceNorthPlumeLoc.png'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
