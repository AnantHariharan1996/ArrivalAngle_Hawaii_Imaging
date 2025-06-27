function vr = variance_reduction( d, df, varargin )
% vr = variance_reduction( d, df, varargin )
%  
% calculated reduction in variance - d is the original data, df is the
% fitted data, varargin is an optional weighting vector.

    if length(varargin) == 0
        % vr = variance_reduction( d, df )
        %
        % d = measured residuals
        % df = from forward problem (df = G*m);

        vr = 100*(1 - (var( d - df ) / var( d )));
    elseif length(varargin) >= 1 
        wt = varargin{1}; 
        vr = 100*(1 - (var( d - df, wt ) / var( d, wt ))); 
    end

end