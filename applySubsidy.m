function [mpc] = applySubsidy(mpc, caseInfo, verbose)
%% applySubsidy: Apply subsidies for certain types of generators
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

    % Set default argin
    if nargin < 3
    	verbose = 1; % show a little debug information
    end

    % Apply subsidies for certain types of generators
    nType = size(caseInfo.subTable{:, 'subType'}, 1);
    for i = 1 : nType
        idx = strcmp(mpc.genfuel, caseInfo.subTable{i, 'subType'});      
        mpc.gencost(idx, 5) = mpc.gencost(idx, 5) - caseInfo.subTable{i, 'subValue'};   
        if verbose == 1
            fprintf('Apply $%d subsidy to gencost for %s\n', ...
                caseInfo.subTable{i, 'subValue'}, caseInfo.subTable{i, 'subType'}{1});
        end    
    end     