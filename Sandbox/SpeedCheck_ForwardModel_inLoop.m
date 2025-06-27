%% Speed check for parallelization purposes
% What does this script do? The goal of
% this is to first benchmark the speed of running the forward model...
clear; clc; close all; 
SyntheticTestFolder =     '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/Stored_ModelSpacePredictions/SyntheticTests/';
HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';
Storagefolder = [HomeDir 'Stored_ModelSpacePredictions/']
% first, declare parameters for model space search. 
Period  = 50;
Storagefolder = [Storagefolder num2str(Period) 's/']; [phvel]=prem_dispersion(1/Period);
cglb = phvel; spacing=0.25; lambda = Period*cglb;
load coastlines  
NoiseFac = 10;
True_Scatterer_Lon = [-160];
True_Scatterer_Lat= [18];
True_Scatterer_Width = [150];
True_Scatterer_Tau = [12];
% 
NewFolder  = [HomeDir 'Raw_ArrivalAngleDeviations/' num2str(Period) 's/'];
Event_Files=  dir([NewFolder '*_elon_elat_lon_lat_phidev_phigc']);
EventNum=1;

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

 LatList = linspace(16,20,100);
ModelCounter=0;
parpool;
tic
plotme = 1;
parfor Modelnum = 1:length(LatList)

     True_Scatterer_Lat = LatList(Modelnum);

% Generate the forward model. 

            % Running forward model! %
        [delta,xgrid,ygrid,ttime_field_perturbed,tau,ttime_field_noscatter,angle,angle_nodiff] ...
            = Get_Arrival_Angle_Residual_GaussianBeam(Period, EvLat,...
            EvLon, True_Scatterer_Lat,True_Scatterer_Lon,True_Scatterer_Tau,True_Scatterer_Width,...
            Ref_YGrid(:),Ref_XGrid(:),cglb,spacing);

        ModelSpaceSearch_Store(Modelnum,:) = delta(:);

%         figure(100)
%         if plotme & length(LatList) < 7
%             subplot(2,3,Modelnum)
% scatter(xgrid,ygrid,50,delta,'filled')
% 
% 
% 
%         end


 end

 toc


 delete(gcp('nocreate'))

 % Elapsed time is 30 seconds for NON-PARFOR execution for 100 iterations
  % Elapsed time is 10 seconds for PARFOR execution