%% Optional: Step 3; compare Different AA
% What does this code do? 
% It uses Foster's mini-array beamforming method to calculate the AAs. 
% Store the AAs and generate a comparison between the
% gradient-based measurements and the beamforming measurements


clear; clc; close all; 

HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';
PredictionsDir =  [HomeDir 'Stored_ModelSpacePredictions/'];
ObservationsDir = [HomeDir 'Raw_ArrivalAngleDeviations/'];
SummaryMisfitDir = [PredictionsDir 'SummaryMisfitStore/'];
SummaryFigDir = [PredictionsDir 'SummaryFigures/'];

EvWeightingFolder = ['/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/Stored_ModelSpacePredictions/EventWeightingFigs/'];
addpath(genpath(pwd))
Min_N_Stations = 3;
L_TOL = 2; % UNITS OF DEGREES DUMBO
Periodlist = [50 66.6667 80 100];
Cmap2Use =     '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/UsefulFunctions/roma.cpt';
load coastlines
Pcounter=0;
Ref_PhVel =4;
for Period = Periodlist
Pcounter=Pcounter+1;
GoodIDFname = [PredictionsDir 'IDLIST_' num2str(Period) 's.txt'];
GoodIDs = load(GoodIDFname);
numevts = length(GoodIDs);

    % Make a massive subplot list 
junk2 = figure(149)
nrows = ceil(sqrt(numevts));
ncols =nrows-1;
if ncols*nrows < numevts
ncols = nrows;
end
ax = subplot_custom_make(junk2,nrows,ncols,[0.05],[0.05],[0.05 0.95],[0.05 0.95]);


    for evnum = 1:length(GoodIDs)
    
        EVID = GoodIDs(evnum);
        AngObsFname = [ObservationsDir num2str(Period) 's/' num2str(EVID)  num2str(Period) 's_elon_elat_lon_lat_phidev_phigc'];
        TTObsFname = [ObservationsDir num2str(Period) 's/' num2str(EVID)  num2str(Period) 's_slon_slat_stt'];
        TTObs_WithAAFname = [ObservationsDir num2str(Period) 's/' num2str(EVID)  num2str(Period) 's_slon_slat_stt_phidev'];

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

        % Now run Foster's approach on the TTime measurements. 
    % 
        [ArrivalAngleList,Best_LocalPhVelList] = ...
    GetArrivalAngles_Event_MiniArray(Slon,Slat,Stt,L_TOL,Elat,Elon,Min_N_Stations);


Ref_X_Grid = [min(Slon) - 2:0.25:max(Slon) + 2];
Ref_Y_Grid = [min(Slat) - 2:0.25:max(Slat) + 2];
[XGRD,YGRD] = meshgrid(Ref_X_Grid,Ref_Y_Grid);

dist2grd_km = deg2km(distance(Elat,Elon,YGRD,XGRD));
TT_Pred  = dist2grd_km./Ref_PhVel;
[ fx,fy,angle,xgrid,ygrid,tgrid2 ] = Get_arrival_angle( Elat,Elon,YGRD(:),XGRD(:),TT_Pred(:),0.25 );

Foster_Angle_Resid = ArrivalAngleList' - griddata(xgrid,ygrid,angle,Slon,Slat);
        % Now Make plots comparing the two approaches
scatter(ax(evnum),Obs_Xgrid,Obs_Ygrid,20,AngleResid,'filled')
hold(ax(evnum),"on")

scatter(ax(evnum),Slon,Slat,100,Foster_Angle_Resid,'filled','MarkerEdgeColor','k','LineWidth',2)

hold(ax(evnum),"on")
plot(ax(evnum),coastlon,coastlat,'color','k','LineWidth',2)
text(ax(evnum),-158, 23,['Event ' num2str(EVID)],'fontsize',16,'FontWeight','bold')

set(ax(evnum),'linewidth',2,'box','on','XTickLabel',[],'YTickLabel',[])
   ax(evnum).XLim=[-159.5 -152.5];
     ax(evnum).YLim=[17 25];
ax(evnum).CLim = [-5 5];
colormap(ax(evnum),jet(20))
 cptcmap(Cmap2Use,ax(evnum),'ncol',20);


TTObsInfo(:,4) = Foster_Angle_Resid;
dlmwrite(TTObs_WithAAFname,TTObsInfo,'delimiter', '\t','precision','%.6f')
clear TTObsInfo
    end

    clear Elonlist
    clear Elatlist


    for eventnum = evnum+1:1:nrows*ncols

ax(eventnum).Color = 'none'
ax(eventnum).Box = 'off'
ax(eventnum).XTick = [];
ax(eventnum).YTick = [];
ax(eventnum).XColor = 'none'
ax(eventnum).YColor = 'none'
ax(eventnum).FontSize = 18

if eventnum == evnum+1
% add colorbar
barbar = colorbar(ax(eventnum),'location','west')
 cptcmap(Cmap2Use,ax(eventnum),'ncol',20);

ax(eventnum).CLim = [-5 5];
ylabel(barbar,{'Arrival Angle','Deviation (˚)'}, 'fontsize',18)
end
    end
text(ax(eventnum),0.5,0.5,[num2str(Period) ' s'],'fontsize',24,'fontweight','bold')

set(junk2,'position',[1637 -77 1134 956])
saveas(junk2,[SummaryFigDir 'OurObsFosterOverlain' num2str(Period) 's' '.png'])
 close all
end