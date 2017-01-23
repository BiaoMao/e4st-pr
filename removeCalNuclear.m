function [mpc, offer] = removeCalNuclear(mpc, offer, caseInfo, verbose)
%% removeCalNuclear: remove nuclear in Califonia
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

    % Extract the bus info
    genBus = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
    genBus = join(genBus, caseInfo.locationInfo);

    idx = strcmp(genBus{:, 'State'}, 'california') & strcmp(mpc.genfuel, 'nuclear');
    [mpc, offer] = removeGen(mpc, offer, idx);

    if verbose == 1
        fprintf('Remove %d nuclear in Califonia\n', length(find(idx)));
    end
end