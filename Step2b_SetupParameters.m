%%% Setup_Parameters

% NOTE,  ONLY FOR RAYLEIGH WAVES!!
[phvel]=prem_dispersion(1/Period)
cglb = phvel;
spacing=0.25;
lambda = Period*cglb;

% Parameter 1: The Position
XsearchList = [-170:1:-154];
YsearchList = [16:1:26];
[XSEARCHGRD,YSEARCHGRD] = meshgrid(XsearchList,YsearchList);
XSEARCHGRD_List=XSEARCHGRD(:); 
YSEARCHGRD_List = YSEARCHGRD(:);

% Parameter 2: The Time Lag
TauMax_List = [40:-2.5:0];

% Parameter 3: The Width of the Anomaly
Loverlambda = [0.5:0.1:2];
L_List = Loverlambda*lambda;


% Implicit Parameter; grid onto which we calculate the model Predictions
%% Actually, this is specified for each event!
% Ref_XGrid = -165:spacing:-145;
% Ref_YGrid = 15:spacing:26;
% [Ref_XGrid,Ref_YGrid] = meshgrid(Ref_XGrid,Ref_YGrid);

N_MODELS = length(TauMax_List)*length(Loverlambda)*length(YSEARCHGRD_List);
    % List of model parameters
    [Ngrid_X,Ngrid_Y,Ngrid_Tau,Ngrid_L] = ndgrid(XsearchList,YsearchList,TauMax_List,L_List);
    Ngrid_X=Ngrid_X(:); Ngrid_Y=Ngrid_Y(:); Ngrid_Tau=Ngrid_Tau(:); Ngrid_L=Ngrid_L(:);
    Lonstore = Ngrid_X;
    Latstore =Ngrid_Y;
    Widthstore= Ngrid_L;
    Taustore=Ngrid_Tau;