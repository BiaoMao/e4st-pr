function [genAuxdata] = updateDamage(busId, genAuxdata, caseInfo, verbose)
%% updateDamage: update locational damages for NOx and SO2
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

    % Set default argin
    if nargin < 4
    	verbose = 1; % show a little debug information
    end

    damInfo = caseInfo.damInfo;
    
    for i = 1 : size(busId, 1)
        idx = damInfo{:, 'bus'} == busId(i);
        damRate = damInfo{idx, 2:3} .* genAuxdata(i, 2:3);
        genAuxdata(i, 5:6) = damRate;
    end

    if verbose == 1
        fprintf('Apply locational damages\n');
    end     
end