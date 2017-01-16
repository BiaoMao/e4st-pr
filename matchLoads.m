function mpc = matchLoads(mpc, caseInfo, verbose)
%% matchLoads: match loads to the real values
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

    % Select full info
    if isfield(caseInfo, 'fullLocInfo')
        locInfo = caseInfo.fullLocInfo;
    else 
        locInfo = caseInfo.locationInfo;
    end

    % Update gen loc infomation
    busTable = array2table(mpc.bus(:, 1), 'VariableNames', {'bus'});
    busTable = join(busTable, locInfo);
    % states = unique(genBus{:, 'State'});
    % statesTable = [table(states, 'VariableNames', {'State'})...
    %     array2table([1 : length(states)]', 'VariableNames', {'load_zone'})];   
    % Sort the State ID first
    caseInfo.loadValue = sortrows(caseInfo.loadValue,'StateID','ascend');
    busTable = join(busTable, caseInfo.loadValue(:, {'State', 'StateID'}));

    % Scale State loads to real values
    load_zone = busTable{:, 'StateID'};    
    opt = struct('pq', 'P', 'scale', 'QUANTITY');
    idxReal = caseInfo.loadValue{:, 'StateID'} ~= 0;
    act_load = caseInfo.loadValue{idxReal, 'realLoads'};
    [mpc.bus, mpc.gen] = scale_load(act_load, mpc.bus, mpc.gen, load_zone,opt);

    % Debug information
    if verbose == 1
        fprintf('Real values of loads are matched\n');
    end
end