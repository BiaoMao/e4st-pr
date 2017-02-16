function mpc = addCapCons(mpc, caseInfo, group, verbose)
%% addCapCons: add capacity constraints for each States and fuel types
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
    
    define_constants;
    capMax = caseInfo.capMax;
    capMin = caseInfo.capMin;
    [m, n] = size(capMax);
    % numCons = m * n;
    numGens = size(mpc.gen, 1);
    %map = zeros(numCons, numGens);
    %cap = zeros(numCons, 1);

    % Select full info
    if isfield(caseInfo, 'fullLocInfo')
        locInfo = caseInfo.fullLocInfo;
    else 
        locInfo = caseInfo.locationInfo;
    end
    % Extract the bus info
    genBus = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
    genBus = join(genBus, locInfo);

    % Select the group to apply these caps
    if strcmp(group, 'new')
        idxGroup = mpc.newgen == 1;
    elseif strcmp(group, 'all')
        idxGroup = ones(numGens, 1);
    end    

    % Set capacity constraints
    idx = 1;
    for i = 1 : m
        for j = 1 : n
            mapIdx = strcmp(genBus{:, 'State'}, capMax.Properties.RowNames{i}) ...
                & strcmp(mpc.genfuel, capMax.Properties.VariableNames{j}) & idxGroup;
            % Check if there is constraint in this group; -1 means no constraints
            if capMax{i, j} == -1
            	continue;
            end
            % Check if there is any gen in this group
            if ~any(mapIdx)
                if capMax{i, j} ~= 0
                    fprintf('No buildable %s exists in %s\n', capMax.Properties.VariableNames{j},...
                        capMax.Properties.RowNames{i});
                end
                continue;
            end
            map(idx, :) = mapIdx;
            % Choose the min value for the targeted value and buildable cap
            if sum(mpc.gen(mapIdx, PMAX)) < capMax{i, j}
                fprintf('Not enough buildable %s exists in %s\n', capMax.Properties.VariableNames{j},...
                     capMax.Properties.RowNames{i});
            end
            caplimMax(idx, 1) = min(capMax{i, j}, sum(mpc.gen(mapIdx, PMAX)));
            caplimMin(idx, 1) = capMin{i, j};
            idx = idx + 1;
        end
    end

    % Update in mpc to make equality constraints
    mpc.caplim.map = map;
    mpc.caplim.max = caplimMax;
    mpc.caplim.min = caplimMin;

    % Debug information
    if verbose == 1
        fprintf('%d Capacity constraints are added\n', size(mpc.caplim.max, 1));
    end
end