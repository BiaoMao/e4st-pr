function [offer] = updateFixedcost(mpc, offer, fuelType, scaling, verbose)
%% updateFixedcost: scale fixed cost by a scaling for new generators
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

    % Adjust cost to build for new generators
    %  Need to change it later, offer(:, 2) = cost2keep + cost2build
    FIXED_COST_IDX = 2;
    idxNew = mpc.newgen == 1;
    nType = length(fuelType);
    for i = 1 : nType
        idx = strcmp(mpc.genfuel, fuelType) & idxNew;      
        offer(idx, FIXED_COST_IDX) = offer(idx, FIXED_COST_IDX) * scaling(i);
        if verbose == 1
            fprintf('Apply subsidy to cost to build for %s\n', caseInfo.subTable{i, 'subType'});
        end    
    end     