% Step 3 Optional

% Record characteristics of the different events, and explore their use for
% weighting. 

clear; clc; close all
HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';
PredictionsDir =  [HomeDir 'Stored_ModelSpacePredictions/'];
ObservationsDir = [HomeDir 'Raw_ArrivalAngleDeviations/'];
SummaryMisfitDir = [PredictionsDir 'SummaryMisfitStore/'];

PeriodList  = [50 66.6667 80 100];

for Period = PeriodList

GoodIDFname = [PredictionsDir 'IDLIST_' num2str(Period) 's.txt'];
GoodIDs = load(GoodIDFname);

    
    for evnum = 1:length(GoodIDs)
    
        EVID = GoodIDs(evnum);
        AngObsFname = [ObservationsDir num2str(Period) 's/' num2str(EVID)  num2str(Period) 's_elon_elat_lon_lat_phidev_phigc'];
        TTObsFname = [ObservationsDir num2str(Period) 's/' num2str(EVID)  num2str(Period) 's_slon_slat_stt'];

        AngObsInfo  = load(AngObsFname,'-ascii');
        TTObsInfo = load(TTObsFname,'-ascii');

        %%%%% Declaring variables
        Elon = AngObsInfo(1,1);
        Elat = AngObsInfo(1,2);
        Obs_Xgrid = AngObsInfo(:,3);
        Obs_Ygrid = AngObsInfo(:,4);
        Elonlist(evnum) = Elon; Elatlist(evnum) = Elat;
        
        AngleResid = AngObsInfo(:,5);
        Weight = AngObsInfo(:,8).*10;
        Slon = TTObsInfo(:,1);
        Slat = TTObsInfo(:,2);
        Stt = TTObsInfo(:,3);
        distkm = deg2km(distance(Elat,Elon,Slat,Slon));
        %%%%%

        % First, make a plot of the locations of the events. 




    end

end