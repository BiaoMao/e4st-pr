function [mpc, offer] = removeNuclear(mpc, offer, busNum, verbose)
%% removeCalNuclear: remove nuclear in busNum
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

    % Remove nuclear in each bus
    for i = 1 : length(busNum)
        curBus = busNum(i);
        idx = strcmp(mpc.genfuel, 'nuclear') & mpc.gen(:, 1) == curBus;
        [mpc, offer] = removeGen(mpc, offer, idx);
    end

    if verbose == 1
        fprintf('Remove %d nuclear\n', length(busNum));
    end
end