function [mpc, offer] = retireGen(mpc, offer, result, caseInfo, group, verbose)
%% retireGen: retire generations which are smaller than theta
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    % Set default argin
    if nargin < 6
        verbose = 1; % show a little debug information
    end

    usedCap = result.reserve.qty.Rp_pos;
            
    % Filter out the small generators
    idx = usedCap < caseInfo.retireTheta;
    usedCap(idx) = 0;

    % Filter out the fuel types that does not retire  
    idxNoretire = zeros(size(mpc.gen, 1), 1);          
    if isfield(caseInfo, 'noRetire')                
        for fuelType = caseInfo.noRetire
            idxNoretire = idxNoretire | strcmp(mpc.genfuel, fuelType);
        end
    end

    % Choose which group to apply retirement
    if strcmp(group, 'new')
        idxNoretire = idxNoretire & (mpc.newgen ~= 1);
    end

    % Set PositiveReserveCap to used capacity
    offer(~idxNoretire, 2) = usedCap;

    % Remove the cost to build,ngt-3,ngcc-9,solar-13,6-wind,, 1-cost to keep    
    isBuilt = mpc.newgen == 1;  
    if ~isempty(mpc.genfuel(isBuilt, :))
        newFuels = unique(mpc.genfuel(isBuilt, :));
        for fuel = newFuels'
            bultIdx = strcmp(mpc.genfuel, fuel);
            offer(bultIdx, 1) = caseInfo.genInfo{fuel, 'Cost2Keep'};
        end
    end

    % Label the gen built in previous decades (newgen - 1)
    mpc.newgen = mpc.newgen - 1;   

    % Delete gen that used caps are zeros    
    mpc.gen(idx, :) = [];
    mpc.gencost(idx, :) = [];
    mpc.genfuel(idx, :) = [];
    mpc.gen_aux_data(idx, :) = [];
    mpc.newgen(idx, :) = [];
    mpc.availability_factor(idx, :) = [];
    offer(idx, :) = [];
    mpc.total_output.map(:, idx) = [];
    mpc.total_output.coeff(idx, :) = [];

    % Debug information
    if verbose == 1
        fprintf('The number of retired generators is %d\n', length(find(idx)));
    end
end