clear; clc; close all; 
% 
%% Step 1: Generate Files Containing 
%% Arrival Angle Measurements
MakePlt = 0;
HomeDir =  '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/';

Periodlist =[50 66.6667 80 100];

for Period = Periodlist
    NewFolder  = [HomeDir 'Raw_ArrivalAngleDeviations/' num2str(Period) 's/'];
    mkdir(NewFolder)

load(['Ma_Data/Ma_Dataset_Output' num2str(Period) 's.mat'])
StanetList = Ma_Dataset_Output.Net;
StanameList = Ma_Dataset_Output.Sta;
Output_Array = Ma_Dataset_Output.Array;
PLUME_DX = find(strcmp(StanetList,'YS'));
Output_Array_PLUME = Output_Array(PLUME_DX,:);
EVENT_IDS = Output_Array_PLUME(:,8);
unique_events = unique(EVENT_IDS);
%% Loop over events. 
for ijk = 1:length(unique_events)

    curr_evID = unique_events(ijk);
%% For every event, output a dtpamp file with the traveltimes. 
    idx = find(EVENT_IDS == curr_evID);
    numstns(ijk) = length(idx);
    if  numstns(ijk) > 15


         current_stalat = Output_Array_PLUME(idx,3);
        current_stalon = Output_Array_PLUME(idx,4);
        current_evlat = Output_Array_PLUME(idx,1);
        current_evlon = Output_Array_PLUME(idx,2);
        current_TTIMES = Output_Array_PLUME(idx,5);
        current_DistanceList = Output_Array_PLUME(idx,11);
        current_EventYear = Output_Array_PLUME(idx,9);
        current_EventDay = Output_Array_PLUME(idx,10);
        titlestr = ['Event: ' num2str(current_EventYear(1)) ', Day:' num2str(current_EventDay(1)) '- Traveltimes'];
        
        
        [ fx,fy,angle,xgrid,ygrid,tgrid2 ] = Get_arrival_angle( current_evlat(1),...
             current_evlon(1),current_stalat,current_stalon,current_TTIMES,0.25 );
        lims = [mean(angle(:))-std(angle(:)) mean(angle(:))+std(angle(:))];
        BestfitC = polyfit(deg2km(current_TTIMES),current_DistanceList,1);
        %% Get arrival angle at locations corresponding to xgrid and ygrid. 
        
        [synth_xlocs] = [min(xgrid(:))-1:0.1:max(xgrid(:))+1];
        [synth_ylocs] = [min(ygrid(:))-1:0.1:max(ygrid(:))+1];
        [SynthXGRD,SynthYGRD] = meshgrid(synth_xlocs,synth_ylocs);
        SynthDists = distance(current_evlat(1),current_evlon(1),SynthYGRD,SynthXGRD);
        
        Pred_TT_SYNTH = deg2km(SynthDists)./BestfitC(1);
        
        [ fx_SYNTH,fy_SYNTH,angle_SYNTH,xgrid_SYNTH,ygrid_SYNTH,tgrid2_SYNTH ] = Get_arrival_angle( current_evlat(1),...
             current_evlon(1),SynthYGRD(:),SynthXGRD(:),Pred_TT_SYNTH(:),0.25 );

Angle_Interped = griddata(xgrid_SYNTH,ygrid_SYNTH,angle_SYNTH,xgrid,ygrid);
Angle_Resid = angle - Angle_Interped;

if MakePlt
        figure(ijk)
subplot(1,2,1)
        scatter(current_stalon,current_stalat,150,current_TTIMES,'filled')
        hold on
        grid on; box on;
        % plot(lontrk,lattrk,'linewidth',2)
        barbar=colorbar;
        ylabel(barbar,'Traveltime (s)')
        colormap(turbo)
title(titlestr)
set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
       load coastlines
       hold on
       plot(coastlon,coastlat,'linewidth',2,'color','k')
       xlim([-162 -148])
       ylim([14 27])
grid on; box on;


subplot(1,2,2)
        scatter(xgrid(:),ygrid(:),50,Angle_Resid(:),'filled')

caxis([-5 5])
   barbar=colorbar;
        ylabel(barbar,'Arrival Angle Deviations (deg)')
        colormap(turbo)
        title('Measured Arrival Angle Residual')
set(gca,'fontsize',18,'fontweight','bold','linewidth',2)
       load coastlines
       hold on
       plot(coastlon,coastlat,'linewidth',2,'color','k')
       xlim([-162 -148])
       ylim([14 27])
grid on; box on;
set(gcf,'position', [266 416 962 395])
saveas(gcf,[NewFolder num2str(curr_evID) num2str(Period) 's' '.png'])
close all
end

[AzCoverageWeight] = Get_Geographical_Weight(current_stalon,...
    current_stalat,xgrid(:),ygrid(:));

%%%%%  Output this- the raw arrival angle 'datum' to be modeled and fit. 
Data2Write(:,1) = current_evlon(1)*ones(size(xgrid(:)));
Data2Write(:,2) = current_evlat(1)*ones(size(xgrid(:)));
Data2Write(:,3) = xgrid(:);
Data2Write(:,4) = ygrid(:);
Data2Write(:,5) = Angle_Resid(:);
Data2Write(:,6) = Angle_Interped(:);
Data2Write(:,7) = tgrid2(:);
Data2Write(:,8) = AzCoverageWeight(:);


RawData2Write(:,1) = current_stalon;
RawData2Write(:,2) =current_stalat;
RawData2Write(:,3) =current_TTIMES;
dlmwrite([NewFolder num2str(curr_evID) num2str(Period) 's' '_elon_elat_lon_lat_phidev_phigc'],Data2Write,'delimiter','\t','precision','%.6f');
dlmwrite([NewFolder num2str(curr_evID) num2str(Period) 's' '_slon_slat_stt'],RawData2Write,'delimiter','\t','precision','%.6f');
clear Data2Write
clear RawData2Write


%%%%%
    end
end
end