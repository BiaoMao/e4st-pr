function mpc = selIsland(mpc, caseInfo, verbose)
%% selIsland: Select islands from MPC to eliminate isolated bus
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

    define_constants;
    island = caseInfo.island;
    % island = 0 means to extract all islands
    if length(island) == 1 && island == 0
    	island = 'all';
    end
    custom.bus{1} = {'bus_island', 'bus_name'};
    custom.gen{1} = {'genfuel', 'gen_aux_data'};
    mpc = extract_islands(mpc, island, custom);

    % Debug information
    if verbose == 1
        if ischar(island)
            fprintf('Extract islands %s from MPC\n', island);
        else
            fprintf('Extract islands %d from MPC\n', island);
        end
    end
end