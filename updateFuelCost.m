function [mpc] = updataFuelCost(mpc, yearInfo, year)
%% updataFuelCost: update fuel cost by the current fuel prices for existing generators
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

	% Format the year as a string
	year = ['Y', num2str(year)];
	idxYear = find(strcmp(yearInfo.Properties.VariableNames, year));

	% For coal
	idx = strcmp(mpc.genfuel, 'coal');
	if isfield(mpc, 'newgen')
		idx = idx & newgen ~= 1; % for old gen
	end
	delta = yearInfo{'coalPrice', idxYear} - yearInfo{'coalPrice', idxYear - 1};
    mpc.gencost(idx, 5) = mpc.gencost(idx, 5) + ...
        delta * mpc.gen_aux_data(idx, 8); % 8th column is the heat rate

    % For oil
	idx = strcmp(mpc.genfuel, 'oil');
	if isfield(mpc, 'newgen')
		idx = idx & newgen ~= 1; % for old gen
	end
	delta = yearInfo{'oilPrice', idxYear} - yearInfo{'oilPrice', idxYear - 1};
    mpc.gencost(idx, 5) = mpc.gencost(idx, 5) + ...
        delta * mpc.gen_aux_data(idx, 8); % 8th column is the heat rate


    % For ng
	idx = strcmp(mpc.genfuel, 'ng') & strcmp(mpc.genfuel, 'ngcc') & strcmp(mpc.genfuel, 'ngt');
	if isfield(mpc, 'newgen')
		idx = idx & newgen ~= 1; % for old gen
	end
	delta = yearInfo{'ngPrice', idxYear} - yearInfo{'ngPrice', idxYear - 1};
    mpc.gencost(idx, 5) = mpc.gencost(idx, 5) + ...
        delta * mpc.gen_aux_data(idx, 8); % 8th column is the heat rate
