%% Step 4: Plot Misfit Surfaces
% What does this script do?
% For every period, we stack the misfit for every event.
% Also, for the lowest misfits for a couple of events, 
% Make plots showing best-fitting model, observations, and 
% Misfit Surface. 
% Also, make plots showing the misfit surfaces for all the reliable events.


clear; clc; close all; load coastlines; 
HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';
PredictionsDir =  [HomeDir 'Stored_ModelSpacePredictions/'];
ObservationsDir = [HomeDir 'Raw_ArrivalAngleDeviations/'];
SummaryMisfitDir = [PredictionsDir 'SummaryMisfitStore/'];
Periodlist = 66.6667;

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


for eventnum = 1:length(RMSList)

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


    % Generate Forward Model for Best-Fit Diffractor, but over a wider
    % space
WiderSpaceXGrid = [-162:0.2:-148];
WiderSpaceYGrid = [13:0.2:27];
[XGRID_WIDE,YGRID_WIDE] = meshgrid(WiderSpaceXGrid,WiderSpaceYGrid);

                % Running forward model! %
        [delta,xgrid,ygrid,ttime_field_perturbed,tau,ttime_field_noscatter,angle,angle_nodiff] ...
            = Get_Arrival_Angle_Residual_GaussianBeam(Period, Elat,...
            Elon, BestLat,BestLon,BestTau,BestWidth,...
            YGRID_WIDE(:),XGRID_WIDE(:),cglb,spacing);
 deltainterp =griddata(xgrid,ygrid,delta,XGRID_WIDE, YGRID_WIDE)      ;
  Step4b_Plotting
    

end


end

