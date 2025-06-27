function ax = subplot_custom_make(varargin)
% ax = subplot_custom_make([figh],Nx,Ny,[fwx],[fwy],[xmami],[ymami])
%
% ARGUMENTS:
% figh  = handle to figure window into which to plot axes [optional]
% Nx    = number of plots in x direction (N columns)
% Ny    = number of plots in y direction (N rows)
% fx    = fraction of x width of subplots to be whitespace between axes [optional]
% fy    = fraction of y height of subplots to be whitespace between axes [optional]
% xmami = [min,max] x extent of figure window within which to fit all axes [optional]
% ymami = [min,max] y extent of figure window within which to fit all axes [optional]

% figure handle in which to plot the new axes
if ishandle(varargin{1})
    figh = varargin{1};
    clf(figh);
    varargin = varargin(2:end);
else
    figh = figure;clf;
end

% set number of plots in each dimension
Nx = varargin{1}; 
Ny = varargin{2};

% set the fraction of x,y extents to be whitespace between axes
if length(varargin)>2
    fwx = varargin{3};
else
    fwx = 0.04*Nx;
end
if length(varargin)>3
    fwy = varargin{4};
else
    fwy = 0.1*Ny;
end

% set the space within the figure window spanned by the axes
if length(varargin)>4
    xmami = varargin{5};
else
    xmami = [0.05 0.99];
end
if length(varargin)>5
    ymami = varargin{6};
else
    ymami = [0.06 0.98];
end



%% okay now make the width/height/spacing_x,spacing_y
hx = diff(xmami)./(Nx + (Nx-1)*fwx);
hy = diff(ymami)./(Ny + (Ny-1)*fwy);
dx = hx*fwx;
dy = hy*fwy; 

%% loop through and make the axes.
% note, y-position is (1) at the top and (end) at the bottom.
%       x-position is (1) on the left and (end) on the right.
for ix = 1:Nx
for iy = 1:Ny
    % ax(Ny+1-iy,Nx+1-ix) = axes(figh,'position',[xmami(1) + (ix-1)*(hx+dx),...
    %                                    ymami(1) + (iy-1)*(hy+dy),...
    %                                    hx,hy],'box','on');
    ax(Ny+1-iy,ix) = axes(figh,'position',[xmami(1) + (ix-1)*(hx+dx),...
                                       ymami(1) + (iy-1)*(hy+dy),...
                                       hx,hy],'box','on');
end
end


end