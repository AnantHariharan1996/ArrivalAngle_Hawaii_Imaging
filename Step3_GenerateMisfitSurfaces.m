%% 
clear; clc; close all;

% For every event, generate misfit surface as a function of variables of interest
HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';
PredictionsDir =  [HomeDir 'Stored_ModelSpacePredictions/'];
ObservationsDir = [HomeDir 'Raw_ArrivalAngleDeviations/'];
SummaryMisfitDir = [PredictionsDir 'SummaryMisfitStore/'];
mkdir(SummaryMisfitDir)
Periodlist =[50] %%% 66.6667 80 100];

Pcounter = 0;
for Period = Periodlist
    Pcounter=Pcounter+1;

    % Now find all the .mat files corresponding to this Period
    FolderName = [PredictionsDir num2str(Period) 's/' ];
    Predictions_List = dir([FolderName 'ModelSpaceSearch_Store_EVID*'])

    Step2b_SetupParameters
    StoreAllMisfits_L1 = zeros(N_MODELS,length(Predictions_List));
    StoreAllMisfits_L2 = zeros(N_MODELS,length(Predictions_List));
    PhVelList = zeros(1,length(Predictions_List));
    RMSList = zeros(1,length(Predictions_List));
    StoreAllVarReduc = zeros(N_MODELS,length(Predictions_List));
    StoreAllMisfits_L1_Weighted = zeros(N_MODELS,length(Predictions_List));
    StoreAllMisfits_L2_Weighted = zeros(N_MODELS,length(Predictions_List));

    for PredCounter = 1:length(Predictions_List)
        disp(['Loaded  ' num2str(100*PredCounter/length(Predictions_List)) '%'])
        currfname = Predictions_List(PredCounter).name;
        full_pred_name = [FolderName currfname];
        load(full_pred_name)

        % For every Model Prediction, get the Event ID and load the
        % observations. 
        
        EVID = [extractBetween(currfname,'EVID','.mat')];
        EVID=EVID{1};

        AngObsFname = [ObservationsDir num2str(Period) 's/' EVID  num2str(Period) 's_elon_elat_lon_lat_phidev_phigc'];
        TTObsFname = [ObservationsDir num2str(Period) 's/' EVID  num2str(Period) 's_slon_slat_stt'];

        AngObsInfo  = load(AngObsFname,'-ascii');
        TTObsInfo = load(TTObsFname,'-ascii');

        %%%%% Declaring variables
        Elon = AngObsInfo(1,1);
        Elat = AngObsInfo(1,2);
        Obs_Xgrid = AngObsInfo(:,3);
        Obs_Ygrid = AngObsInfo(:,4);

        
        AngleResid = AngObsInfo(:,5);
        Weight = AngObsInfo(:,8).*10;
        Slon = TTObsInfo(:,1);
        Slat = TTObsInfo(:,2);
        Stt = TTObsInfo(:,3);
        distkm = deg2km(distance(Elat,Elon,Slat,Slon));
        %%%%%

        % get best-fit phase velocity, assuming gc propagation
         [p,S] = polyfit(distkm,Stt,1); gcphvel = 1/p(1);
        % predict ttimes
        ttpred =polyval(p,distkm);
        % get tt error
        errors = abs(ttpred-Stt); avgerror = mean(errors);

         % Now, loop over model predictions and calculate misfit for this
         % event. 
          Lonstore = Output_Store.Lonstore;
          Latstore = Output_Store.Latstore;
          Widthstore = Output_Store.Widthstore;
          Taustore = Output_Store.Taustore;



          
           xgrid =Output_Store.xgrid;
           ygrid =Output_Store.ygrid;
           
           L1_MisfitStore = 999999.*ones(size(Taustore)); 
           L2_MisfitStore = L1_MisfitStore;
           
           L1_MisfitStore_Weighted = L2_MisfitStore;
           L2_MisfitStore_Weighted= L2_MisfitStore;

           Var_Reduc_List = L1_MisfitStore;
          for modelspace_num = 1:length(Taustore)
              current_predictions = Output_Store.ModelSpaceSearch_Store(modelspace_num,:);

              %  interpolate the predictions on the observations' locations. 

              Interped_AA = griddata(xgrid,ygrid,current_predictions,Obs_Xgrid,Obs_Ygrid);
               
                Difference = AngleResid-Interped_AA;
                Mean_L1Misfit = mean(abs(Difference)); % currently unweighted (but still avg) misfit
                Mean_L2Misfit = mean(abs(Difference.^2)); % currently unweighted (but still avg) misfit
                % Now do weighted calculations
                AbsDiff = abs(Difference); absDiffsquared = abs(Difference.^2);
                WmeanL1  = sum(AbsDiff.*Weight)./sum(Weight);
                WmeanL2  = sum(absDiffsquared.*Weight)./sum(Weight);

                L1_MisfitStore(modelspace_num) = Mean_L1Misfit;
                L2_MisfitStore(modelspace_num) = Mean_L2Misfit;
                L1_MisfitStore_Weighted(modelspace_num) = WmeanL1;
                L2_MisfitStore_Weighted(modelspace_num)= WmeanL2;
          
                % Get Variance Reduction
                Var_Reduc_List(modelspace_num) =  variance_reduction( AngleResid',Interped_AA');
          end
       
         StoreAllMisfits_L1(:,PredCounter) = L1_MisfitStore;
         StoreAllMisfits_L2(:,PredCounter) = L2_MisfitStore;
         StoreAllVarReduc(:,PredCounter) = Var_Reduc_List;
         PhVelList(PredCounter) = gcphvel;
        RMSList(PredCounter) = avgerror;
        EVIDLIST{PredCounter} = EVID;
       
         StoreAllMisfits_L1_Weighted(:,PredCounter) = L1_MisfitStore_Weighted;
         StoreAllMisfits_L2_Weighted(:,PredCounter) = L2_MisfitStore_Weighted;



    end


         Lonstore = Output_Store.Lonstore;
        Latstore = Output_Store.Latstore;
        Widthstore = Output_Store.Widthstore;
        Taustore = Output_Store.Taustore;



    MisfitSurfaceSummary.StoreAllMisfits_L1=StoreAllMisfits_L1;
    MisfitSurfaceSummary.StoreAllMisfits_L2=StoreAllMisfits_L2;
    MisfitSurfaceSummary.PhVelList=PhVelList;
    MisfitSurfaceSummary.RMSList=RMSList;
    MisfitSurfaceSummary.EVIDLIST=EVIDLIST;
    MisfitSurfaceSummary.Var_Reduc_List=StoreAllVarReduc;
    MisfitSurfaceSummary.Lonstore=Lonstore;
    MisfitSurfaceSummary.Latstore=Latstore;
    MisfitSurfaceSummary.Widthstore=Widthstore;
    MisfitSurfaceSummary.Taustore=Taustore;
    MisfitSurfaceSummary.StoreAllMisfits_L1_Weighted=StoreAllMisfits_L1_Weighted;
    MisfitSurfaceSummary.StoreAllMisfits_L2_Weighted=StoreAllMisfits_L2_Weighted;
    save([SummaryMisfitDir '/MisfitSurfaces_' num2str(Period) 's.mat'],'MisfitSurfaceSummary') 
clear MisfitSurfaceSummary
clear StoreAllMisfits_L1; 
clear StoreAllMisfits_L2
clear PhVelList
clear RMSList
clear EVIDLIST
clear Var_Reduc_List
clear StoreAllMisfits_L1_Weighted; 
clear StoreAllMisfits_L2_Weighted
end

