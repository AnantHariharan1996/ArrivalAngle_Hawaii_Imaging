function [ fx,fy,angle,xgrid,ygrid,tgrid2 ] = Get_arrival_angle( evla,evlo,stlat,stlon,tsec,spacing )
eqlat=evla(1);
eqlon=evlo(1);
%d,az]=distance(eqlat,eqlon,stlat,stlon);
% 
% dkm=deg2km(distance(eqlat,eqlon,stlat,stlon));
% 
% p=polyfit(dkm,tsec,1);
% ave_c=1./p(1)
% tp=polyval(p,dkm);
% meanv=mean(abs(tp-tsec));
% h=find(abs(tp-tsec)>(3*meanv)); %remove outlier travel times
% stlat(h)=[];
% stlon(h)=[];
% tsec(h)=[];

xmin=min(stlon);
xmax=max(stlon);
ymin=min(stlat);
ymax=max(stlat);
[xgrid,ygrid] = meshgrid(xmin:spacing:xmax,ymin:spacing:ymax);
sincolat = sind(90-ygrid);

%interpolate travel times and amplitudes onto grid
tgrid=griddata(stlon,stlat,tsec,xgrid,ygrid);
newt=tgrid(:);

tgrid2=reshape(newt,[size(xgrid)]);

[fx,fy]=gradient(tgrid2,deg2km(spacing),deg2km(spacing)); %dt/dtheta and dt/dphi

angle = rad2deg(atan2((fx./sincolat),fy));
%angle = rad2deg(atan2((fx),fy));

idx = find(angle < 0);
angle(idx) = 360+angle(idx);
  
k = convhull(stlon,stlat);
polylon=stlon(k);
polylat=stlat(k);

in = inpolygon(xgrid,ygrid,polylon,polylat);
fx=fx(in);
fy=fy(in);
angle=angle(in);
xgrid=xgrid(in);
ygrid=ygrid(in);
tgrid2=tgrid2(in);

nonan=find(isnan(angle)==0);
fx=fx(nonan);
fy=fy(nonan);
angle=angle(nonan);
xgrid=xgrid(nonan);
ygrid=ygrid(nonan);
tgrid2=tgrid2(nonan);
end

