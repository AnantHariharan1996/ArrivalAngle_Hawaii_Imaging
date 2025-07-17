M = 50; %number of independent observations; number of fds if you are only fitting fds
N = 3; % number of model parameters; columns in Ftest table (numerator degrees of freedom, df1)
% M-N corresponds to the rows in the Ftest table (denominator degrees of freedom, df2)
%in95 =1; % Set outside this script...
FDT = load('F_Distribution_Tables.mat'); 

in95 = 0
if in95==1
    Confidence_Level    = 95;
else
    Confidence_Level    = 90;
end

Significance_Level  = num2str(100-Confidence_Level);
F_Dist_Table    = FDT.F_Distribution_Tables.(['Alpha_' Significance_Level]);
ftest = (1+((N/(M-N))*F_Dist_Table(M-N,N)))^0.5 %%%