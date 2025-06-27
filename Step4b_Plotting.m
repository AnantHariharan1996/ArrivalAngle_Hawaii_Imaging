Cmap2Use= '/Users/ananthariharan/Documents/GitHub/ArrivalAngle_Hawaii_Imaging/UsefulFunctions/roma.cpt';  
%%%%%%%%%%%% Make Figures %%%%%%%%%%%
    figure(eventnum)
    % First plot observations
    ax1=subplot(2,4,[1 2])
    scatter(Obs_Xgrid,Obs_Ygrid,50,AngleResid,'filled')
hold on
    barbar=colorbar;
    ylabel(barbar,'Arrival Angle Deviation (°)')
    caxis([-5 5])
    colormap(turbo)
        hold on
    plot(coastlon,coastlat,'linewidth',2,'color','k')
    scatter(RawSlon,RawSlat,50,[1 0 1],'^','filled')

    box on
    set(gca,'linewidth',2,'fontsize',17,'layer','top','fontweight','bold')
    text(-160,16,'Observations','fontweight','bold','fontsize',17)
           xlim([-162 -148])
       ylim([14 27])
    % Second, plot the best-fitting predictions
    ax2=subplot(2,4,[3 4])

        contourf(XGRID_WIDE,YGRID_WIDE,deltainterp,50,'EdgeColor','none')
        hold on
            text(-160,16,'Predictions','fontweight','bold','fontsize',17)

               viscircles(ax2,[BestLon BestLat],km2deg(BestWidth),'color','r','linewidth',4)
    caxis([-5 5])
    colormap(turbo)
        hold on
    plot(coastlon,coastlat,'linewidth',2,'color','k')
        barbar=colorbar;
    ylabel(barbar,'Arrival Angle Deviation (°)')
        box on
    set(gca,'linewidth',2,'fontsize',17,'layer','top','fontweight','bold')
           xlim([-162 -148])
       ylim([14 27])
    % Third, plot the Misfit as function of Width
   ax3= subplot(2,4,[5])
    plot(VaryingWidths_BestSection,L1Misfit_VaryingWidths,'-ro','linewidth',2)
    % Fourth, plot the Misfit as function of Tau
        ylabel('L1 Misfit (degrees)')
    xlabel('Width (km)')
        box on
    set(gca,'linewidth',2,'fontsize',17,'layer','top','fontweight','bold')
   ax4= subplot(2,4,[6])
    plot(VaryingTaus_BestSection,L1Misfit_VaryingTaus,'-ro','linewidth',2)
    ylabel('L1 Misfit (degrees)')
    xlabel('\tau (s)')
    box on
    set(gca,'linewidth',2,'fontsize',17,'layer','top','fontweight','bold')
    % Fifth, plot the Misfit as function of Position
   ax5= subplot(2,4,[7 8])
    scatter(VaryingLons_BestSection,VaryingLats_BestSection,200,L1Misfit_VaryingLonLats,'filled','o')
    barbar=colorbar;
    ylabel(barbar,'L1 Misfit (degrees)')
    hold on
    plot(coastlon,coastlat,'linewidth',2,'color','k')
    box on
    set(gca,'linewidth',2,'fontsize',17,'layer','top','fontweight','bold')
       xlim([-164 -153])
       ylim([15 27])
       set(figure(eventnum),'Position', [146 149 958 651])
       sgtitle(titlestr,'fontweight','bold','fontsize',20)
% throw on a different colormap 

 cptcmap(Cmap2Use,ax5,'ncol',20);

set(gcf,'position',[1519 172 1002 636])

saveas(figure(eventnum),FigureName)
close all