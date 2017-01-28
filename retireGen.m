function [mpc, offer, result] = retireGen(mpc, offer, result, caseInfo, group, verbose)
%% retireGen: retire generations which are smaller than theta or not used
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
    
    % Calculate and save shutdown capacity
    usedCap = result.reserve.qty.Rp_pos;
    idxExist = (mpc.newgen ~= 1);
    idxNew = mpc.newgen == 1;
    capShutdown = sum(offer(idxExist, 2) - usedCap(idxExist, 1));
            
    % Filter out the small generators
    idxSmall = usedCap < caseInfo.retireTheta;   
    usedCap(idxSmall) = 0;

    % Extract the bus info
    genBus = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
    genBus = join(genBus, caseInfo.locationInfo);

    % Choose which group to apply retirement
    % Only apply to generators
    if strcmp(group, 'new')
        idxRetire = idxNew;
    elseif strcmp(group, 'ca-or-new')
        idxRetire = idxNew |...
            (strcmp(genBus{:, 'Nation'}, 'CA') & ~strcmp(mpc.genfuel, 'dl'));
    elseif strcmp(group, 'all')
        idxRetire = ~strcmp(mpc.genfuel, 'dl');
    end

    % Filter out the fuel types that does not retire  
    idxNoretire = zeros(size(mpc.gen, 1), 1);          
    if isfield(caseInfo, 'noRetire')                
        for fuelType = caseInfo.noRetire'
            idxNoretire = idxNoretire | strcmp(mpc.genfuel, fuelType);
        end
        % Only apply no retire to existing gens
        idxNoretire = idxNoretire & idxExist; 
    end  
    idxRetire = idxRetire & (~idxNoretire);

    % Retirement : set PositiveReserveCap to used capacity       
    offer(idxRetire, 2) = usedCap(idxRetire, 1);
    mpc.gen(idxRetire, 9) = usedCap(idxRetire, 1); % PMAX

    % Remove the cost to build for new generators  
    isBuilt = mpc.newgen == 1;  
    if ~isempty(mpc.genfuel(isBuilt, :))
        newFuels = unique(mpc.genfuel(isBuilt, :));
        for fuel = newFuels'
            bultIdx = strcmp(mpc.genfuel, fuel) & isBuilt;
            % For oswind: fixed cost is 0 after beening built
            if strcmp(fuel, 'oswind')                
                offer(bultIdx, 1) = 0;
                continue;
            end
            % For other built fuel types
            offer(bultIdx, 1) = offer(bultIdx, 1) - caseInfo.genInfo{fuel, 'Cost2Build'}; 
        end
    end

    % Label the gen built in previous decades (newgen - 1)
    if strcmp(group, 'new')
        mpc.newgen = mpc.newgen - 1; 
    elseif strcmp(group, 'all')            
        mpc.newgen = mpc.newgen - 1; 
    end  

    % Delete generators that caps are zeros    
    % idxSmall includes zero caps
    idx = idxSmall & idxRetire;    
    [mpc, offer] = removeGen(mpc, offer, idx);
    result = removeRes(result, idx);

    % Debug information
    if verbose == 1
        fprintf('The number of retired generators is %d\n', length(find(idx)));
        fprintf('The capacity of shutdown generators is %f\n', capShutdown);
    end
end