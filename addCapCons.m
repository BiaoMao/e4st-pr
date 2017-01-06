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

    capConstraints = caseInfo.capConstraints;
    [m, n] = size(capConstraints);
    numCons = m * n;
    numGens = size(mpc.gen, 1);
    map = zeros(numCons, numGens);
    cap = zeros(numCons, 1);

    % Select the group to apply these caps
    if strcmp(group, 'new')
        idxGroup = mpc.newgen == 1;
    elseif strcmp(group, 'all')
        idxGroup = ones(numGens, 1);
    end    

    % Extract the bus info
    genBus = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
    genBus = join(genBus, caseInfo.locationInfo);
    idx = 1;
    for i = 1 : m
        for j = 1 : n
            map(idx, :) = strcmp(genBus{:, 'State'}, capConstraints.Properties.RowNames{i}) ...
                & strcmp(mpc.genfuel, capConstraints.Properties.VariableNames{j}) & idxGroup;
            cap(idx, 1) = capConstraints{i, j};
            idx = idx + 1;
        end
    end

    % Update in mpc to make equality constraints
    mpc.caplim.map = map;
    mpc.caplim.max = cap;
    %mpc.caplim.min = cap;

    % Debug information
    if verbose == 1
        fprintf('Capacity constraints are added\n');
    end
end