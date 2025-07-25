addpath(genpath(pwd))
clear; clc; close all; 
% 
%% Step 2: Generate Files Containing 
%% Model Space Predictions for Every Event in Dataset. 

HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';

Match_Periodlist =[50 66.6667 80 100];
Match_Event_Search_Limits = [12189 10652 7357 7795];

PCounter = 0;
Periodlist = [50]

for Period = Periodlist
    
    PCounter = PCounter +1;
    current_eventlimit = Match_Event_Search_Limits(find(Match_Periodlist== Period));
    NewFolder  = [HomeDir 'Raw_ArrivalAngleDeviations/' num2str(Period) 's/'];
    NewFolder_ModelPredictions  = [HomeDir 'Stored_ModelSpacePredictions/' num2str(Period) 's/'];
    NewFolder_ModelFigs  = [HomeDir 'Stored_ModelSpacePredictions/' num2str(Period) 's_Figs/'];
    mkdir(NewFolder_ModelPredictions); mkdir(NewFolder_ModelFigs);
    Step2b_SetupParameters

     % Loop over Files
     Event_Files=  dir([NewFolder '*_elon_elat_lon_lat_phidev_phigc']);

     for ijk = 1:1:length(Event_Files) %:-1:1

           Current_Event_File = Event_Files(ijk).name;
           % Get Event ID and geometry
           CurrID =extractBefore(Current_Event_File,[num2str(Period) 's']);
           CurrID_num =  str2num(CurrID); CurrInfo = load([NewFolder Current_Event_File],'-ascii');
           EvLon = CurrInfo(1,1); EvLat = CurrInfo(1,2);
           CurrRawDataInfo = load([NewFolder CurrID  num2str(Period) 's' '_slon_slat_stt'],'-ascii');
           FullFname_Current_Event_File = [NewFolder Current_Event_File];
            
           ObsInfo = load(FullFname_Current_Event_File,'-ascii')
           Obs_Lon = ObsInfo(:,3);  Obs_Lat = ObsInfo(:,4); Obs_Delta = ObsInfo(:,5);

         Ref_XGrid = min(Obs_Lon)-0.25:spacing:max(Obs_Lon)+0.25;
         Ref_YGrid = min(Obs_Lat)-0.25:spacing:max(Obs_Lat)+0.25;
         
        % Set up the grid automatically based on the observations. 
        [Ref_XGrid,Ref_YGrid] = meshgrid(Ref_XGrid,Ref_YGrid);

        if Period == 80
            if CurrID_num > 8268 & CurrID_num < 19002
                CurrID_num
           Step2c_RunModelSearch
            else
            disp('This must be an event in phase 2 of the PLUME deployment.')
            end

        else
           % Loop over Model Parameters and store the predictions here. 
            if CurrID_num <= current_eventlimit
                CurrID_num
           Step2c_RunModelSearch
            else
            disp('This must be an event in phase 2 of the PLUME deployment.')
            end
        end

     end
            disp('Done all events for this period.')

end