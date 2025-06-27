function [AzCoverageWeight] = Get_Geographical_Weight(stalon,stalat,xgrid,ygrid)
% Calculate the number of azimuthal bins that are sampled from all
% stations within 1 degrees
AzCoverageWeight = zeros(size(xgrid));
for junkcounter= 1:length(xgrid)

    [dist2stns,az2stns] = distance(ygrid(junkcounter),xgrid(junkcounter),stalat,stalon);
    idx = find(dist2stns < 1);
    closestalons = stalon(idx); closestalats = stalat(idx);
    close_azs = az2stns(idx);
    % Now find how many azimuthal bins this corresponds to
    
    azcounter = 0;
    Bin1Dx = (find(close_azs > 0 & close_azs <= 90));
    Bin2Dx = (find(close_azs > 90 & close_azs <= 180));
    Bin3Dx = (find(close_azs > 180 & close_azs <= 270));
    Bin4Dx = (find(close_azs > 270 & close_azs <= 360));
    
    if length(Bin1Dx) > 0
    azcounter =     azcounter +1;
    end
    if length(Bin2Dx) > 0
    azcounter =     azcounter +1;
    end
    if length(Bin3Dx) > 0
    azcounter =     azcounter +1;
    end
    if length(Bin4Dx) > 0
    azcounter =     azcounter +1;
    end

    AzCoverageWeight(junkcounter) = azcounter;
end



end