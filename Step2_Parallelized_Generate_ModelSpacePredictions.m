addpath(genpath(pwd))
clear; clc; close all; 
% 
%% Step 2: Generate Files Containing 
%% Model Space Predictions for Every Event in Dataset. 
load coastlines

HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';

Match_Periodlist =[50 66.6667 80 100];
Match_Event_Search_Limits = [12189 10652 7357 7795];

PCounter = 0;
Periodlist = [50];
delete(gcp('nocreate'))
parpool; 

%% Declaring some variables
MatSuffix = '.mat';
NameStringMiddle = 'ModelSpaceSearch_Store_EVID';
% DoNotOverWrite; if this variable is set to 1 then skips if already done.
DoNotOverWrite = 0;
%%
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

     parfor ijk = 1:1:length(Event_Files) %:-1:1

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
           % Loop over Model Parameters and store the predictions here. 
            if CurrID_num <= current_eventlimit
                CurrID_num
%% Step 2c

% Grid search over Tau, ScattererLocations, and Width. 
OutputFname = [ NewFolder_ModelPredictions NameStringMiddle CurrID MatSuffix];

if exist(OutputFname,'File') == 2 && DoNotOverWrite == 1
    % do nothing
    disp(['Predictions for Event ' CurrID ' already exist!'])
else 
    % execute grid search
maxnumvals = length(TauMax_List)*length(L_List)*length(XSEARCHGRD_List);

%Obs_Lon = ObsInfo(:,3);  Obs_Lat = ObsInfo(:,4); Obs_Delta = ObsInfo(:,5);


ModelSpaceSearch_Store = zeros([maxnumvals length(Ref_XGrid(:))]);

    Lonstore = Ngrid_X;
    Latstore =Ngrid_Y;
    Widthstore= Ngrid_L;
    Taustore=Ngrid_Tau;
  
        for ModelCounter = 1:length(Taustore)
             [ num2str(Period) 's, Event ' CurrID ': '    num2str(100*ModelCounter/maxnumvals) '% Complete']
             
            scattererlon = Lonstore(ModelCounter);
            scattererlat = Latstore(ModelCounter);
current_timelag = Taustore(ModelCounter);
current_width  = Widthstore(ModelCounter);
            % Running forward model! %
        [delta,xgrid,ygrid,ttime_field_perturbed,tau,ttime_field_noscatter,angle,angle_nodiff] ...
            = Get_Arrival_Angle_Residual_GaussianBeam(Period, EvLat,...
            EvLon, scattererlat,scattererlon,current_timelag,current_width,...
            Ref_YGrid(:),Ref_XGrid(:),cglb,spacing);
        


        ModelSpaceSearch_Store(ModelCounter,:) = delta(:);



%% Occasionally, save a figure showing results from the forward model
if rem(ModelCounter,2000) == 1000

    OutputFigname = [ NewFolder_ModelFigs 'EVID_'...
     CurrID '_' num2str(current_timelag) '_tau_' ...
     num2str(scattererlat) '_' num2str(scattererlon) '_latlon_' num2str(current_width) '_width.png'];
    
    figure(ModelCounter+ijk)
    subplot(1,2,1)
    scatter(xgrid(:),ygrid(:),50,delta(:),'filled')
    barbar=colorbar;
    caxis([-5 5])
    hold on
    scatter(scattererlon,scattererlat,15,'k','filled')
    title({['Scatterer Lon: ' num2str(scattererlon) ', Scatterer Lat: ' num2str(scattererlat) ],['\tau = ' num2str(current_timelag),', width: ' num2str(current_width) ' km']})
    hold on
    plot(coastlon,coastlat,'linewidth',2,'color','k')
    box on; 
    colormap(flipud(turbo))
    set(gca,'linewidth',2,'fontsize',16)
    xlim([min(xgrid(:)) max(xgrid(:))])
    ylim([min(ygrid(:)) max(ygrid(:))])
    ylabel(barbar,'Arrival Angle Deviation (degrees)')
  

        subplot(1,2,2)
            scatter(Obs_Lon(:),Obs_Lat(:),50,Obs_Delta(:),'filled')
    barbar=colorbar;
    caxis([-5 5])
    hold on
   plot(coastlon,coastlat,'linewidth',2,'color','k')
    box on; 
    colormap(flipud(turbo))
    set(gca,'linewidth',2,'fontsize',16)
    xlim([min(xgrid(:)) max(xgrid(:))])
    ylim([min(ygrid(:)) max(ygrid(:))])
    ylabel(barbar,'Arrival Angle Deviation (degrees)')
hold on
scatter(CurrRawDataInfo(:,1),CurrRawDataInfo(:,2),150,[1 0 1],'^','filled')
title(['Observed Arrival Angles: EVID ' CurrID])
    set(figure(ModelCounter+ijk),'position',[129 395 791 317])



    saveas(figure(ModelCounter+ijk),OutputFigname)
    % 
end


%%


        end
    
    
    




         Output_Store = struct("ModelSpaceSearch_Store",ModelSpaceSearch_Store,...
             "Lonstore",Lonstore,...
             "Latstore",Latstore,"Widthstore",Widthstore,...
             "Taustore",Taustore,"xgrid",xgrid,"ygrid",ygrid);

parsave_OutputStore(OutputFname,Output_Store)
end
            
            
            
            else
            disp('This must be an event in phase 2 of the PLUME deployment.')
            end
        

     end
            disp('Done all events for this period.')

end