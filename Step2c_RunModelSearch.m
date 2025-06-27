%% Step 2c

% Grid search over Tau, ScattererLocations, and Width. 
OutputFname = [ NewFolder_ModelPredictions 'ModelSpaceSearch_Store_EVID' CurrID '.mat'];

% DoNotOverWrite; if this variable is set to 1 then skips if already done.
DoNotOverWrite = 0;

if exist(OutputFname) == 2 && DoNotOverWrite == 1
    % do nothing
    disp(['Predictions for Event ' CurrID ' already exist!'])
else 
    % execute grid search
ModelCounter=0;
load coastlines
clear ModelSpaceSearch_Store Lonstore Latstore Taustore Widthstore
maxnumvals = length(TauMax_List)*length(L_List)*length(XSEARCHGRD_List);

%Obs_Lon = ObsInfo(:,3);  Obs_Lat = ObsInfo(:,4); Obs_Delta = ObsInfo(:,5);


ModelSpaceSearch_Store = zeros([maxnumvals length(Ref_XGrid(:))]);
for taucount = 1:length(TauMax_List)
    current_timelag = TauMax_List(taucount);

    
    for widthcount = 1:length(L_List)
        current_width = L_List(widthcount);

        for positioncount = 1:length(XSEARCHGRD_List)
             [ num2str(Period) 's, Event ' CurrID ': '    num2str(100*ModelCounter/maxnumvals) '% Complete']
             
            scattererlon = XSEARCHGRD_List(positioncount);
            scattererlat = YSEARCHGRD_List(positioncount);

            % Running forward model! %
        [delta,xgrid,ygrid,ttime_field_perturbed,tau,ttime_field_noscatter,angle,angle_nodiff] ...
            = Get_Arrival_Angle_Residual_GaussianBeam(Period, EvLat,...
            EvLon, scattererlat,scattererlon,current_timelag,current_width,...
            Ref_YGrid(:),Ref_XGrid(:),cglb,spacing);
        

        ModelCounter=ModelCounter+1;

        ModelSpaceSearch_Store(ModelCounter,:) = delta(:);
        Lonstore(ModelCounter) =scattererlon;
        Latstore(ModelCounter) =scattererlat;
        Widthstore(ModelCounter)= current_width;
        Taustore(ModelCounter)=current_timelag;


%% Occasionally, save a figure showing results from the forward model
if rem(ModelCounter,2000) == 1000

    OutputFigname = [ NewFolder_ModelFigs 'EVID_'...
     CurrID '_' num2str(current_timelag) '_tau_' ...
     num2str(positioncount) '_poscount_' num2str(current_width) '_width.png'];
    
    figure(ModelCounter)
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
    set(gcf,'position',[129 395 791 317])



    saveas(gcf,OutputFigname)
    close all
    % 
end


%%


        end
    
    
    end


end

        Output_Store.ModelSpaceSearch_Store = ModelSpaceSearch_Store;
        Output_Store.Lonstore = Lonstore;
        Output_Store.Latstore = Latstore;
        Output_Store.Widthstore = Widthstore;
        Output_Store.Taustore = Taustore;
        Output_Store.xgrid = xgrid(:);
        Output_Store.ygrid = ygrid(:);
        
save([ NewFolder_ModelPredictions 'ModelSpaceSearch_Store_EVID' CurrID '.mat'],'Output_Store')
clear Output_Store
end