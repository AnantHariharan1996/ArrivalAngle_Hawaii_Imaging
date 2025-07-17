%% Store list of 'good' events for each period of interest
%
clear; clc; close all;


HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';

Match_Periodlist =[50 66.6667 80 100];
Match_Event_Search_Limits = [12189 10652 7357 7795];

PCounter = 0;
Periodlist = Match_Periodlist;  % [50]

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
    IDCounter = 0;
    WIDE_IDCounter=0;
     for ijk = 1:1:length(Event_Files) %:-1:1

           Current_Event_File = Event_Files(ijk).name;
           % Get Event ID and geometry
           CurrID =extractBefore(Current_Event_File,[num2str(Period) 's']);
           CurrID_num =  str2num(CurrID); CurrInfo = load([NewFolder Current_Event_File],'-ascii');
           EvLon = CurrInfo(1,1); EvLat = CurrInfo(1,2);
           CurrRawDataInfo = load([NewFolder CurrID  num2str(Period) 's' '_slon_slat_stt'],'-ascii');
           FullFname_Current_Event_File = [NewFolder Current_Event_File];
            
           ObsInfo = load(FullFname_Current_Event_File,'-ascii');
           Obs_Lon = ObsInfo(:,3);  Obs_Lat = ObsInfo(:,4); Obs_Delta = ObsInfo(:,5);

         Ref_XGrid = min(Obs_Lon)-0.25:spacing:max(Obs_Lon)+0.25;
         Ref_YGrid = min(Obs_Lat)-0.25:spacing:max(Obs_Lat)+0.25;
         
        % Set up the grid automatically based on the observations. 
        [Ref_XGrid,Ref_YGrid] = meshgrid(Ref_XGrid,Ref_YGrid);


        % Store IDs in list of events. 

           if Period == 50 && CurrID_num <= 12189
                IDCounter = IDCounter+1;
                IDLIST{IDCounter} = num2str(CurrID_num);

           elseif Period == 66.6667  && CurrID_num <= 10652
                IDCounter = IDCounter+1;
                IDLIST{IDCounter} = num2str(CurrID_num);

           elseif Period == 80 && (CurrID_num <= 7357 || (CurrID_num > 8268 && CurrID_num < 19002))
                IDCounter = IDCounter+1;
                IDLIST{IDCounter} = num2str(CurrID_num);

           elseif Period == 100 && CurrID_num <= 7795
                IDCounter = IDCounter+1;
                IDLIST{IDCounter} = num2str(CurrID_num);


           end


           if Period == 50 && CurrID_num  > 12189
                WIDE_IDCounter = WIDE_IDCounter+1;
                WIDEIDLIST{WIDE_IDCounter} = num2str(CurrID_num);

           elseif Period == 66.6667  && CurrID_num >  10652
                WIDE_IDCounter = WIDE_IDCounter+1;
                WIDEIDLIST{WIDE_IDCounter} = num2str(CurrID_num);

           elseif Period == 80 && (CurrID_num > 7357 || (CurrID_num < 8268 && CurrID_num > 19002))
                WIDE_IDCounter = WIDE_IDCounter+1;
                WIDEIDLIST{WIDE_IDCounter} = num2str(CurrID_num);

           elseif Period == 100 && CurrID_num  >  7795
                WIDE_IDCounter = WIDE_IDCounter+1;
                WIDEIDLIST{WIDE_IDCounter} = num2str(CurrID_num);


           end





     end


    % Write out the list of good events. 
    IDs_fname = [HomeDir 'Stored_ModelSpacePredictions/IDLIST_' num2str(Period) 's.txt' ];
    writecell(IDLIST',IDs_fname)
clear IDLIST


    % Write out the list of good events. 
    WIDE_IDs_fname = [HomeDir 'Stored_ModelSpacePredictions/WIDE_IDLIST_' num2str(Period) 's.txt' ];
    writecell(WIDEIDLIST',WIDE_IDs_fname)
clear WIDEIDLIST



end