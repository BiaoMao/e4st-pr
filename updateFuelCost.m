function [mpc] = updateFuelCost(mpc, yearInfo, year, verbose)
%% updataFuelCost: update fuel cost by the current fuel prices for existing generators
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

	% Format the year as a string
	year = ['Y', num2str(year)];
	idxYear = find(strcmp(yearInfo.Properties.VariableNames, year));

	% For coal: lower limit is VOM
	idx = strcmp(mpc.genfuel, 'coal');
	if isfield(mpc, 'newgen')
		idx = idx & mpc.newgen ~= 1; % for non-new generators
	end
	delta = yearInfo{'coalPrice', idxYear} - yearInfo{'coalPrice', idxYear - 1};	
    mpc.gencost(idx, 5) = max(mpc.gencost(idx, 5) + ...
        delta * mpc.gen_aux_data(idx, 8), mpc.gen_aux_data(idx, 9)); % 8th column is the heat rate, 9th is VOM
    if verbose == 1
    	fprintf('The coal price has changed %.2f in %s for %d gens\n', delta, year, length(find(idx)))
    end        

    % For oil
	idx = strcmp(mpc.genfuel, 'oil');
	if isfield(mpc, 'newgen')
		idx = idx & mpc.newgen ~= 1; % for non-new generator
	end
	delta = yearInfo{'oilPrice', idxYear} - yearInfo{'oilPrice', idxYear - 1};
    mpc.gencost(idx, 5) = max(mpc.gencost(idx, 5) + ...
        delta * mpc.gen_aux_data(idx, 8), mpc.gen_aux_data(idx, 9)); % 8th column is the heat rate
    if verbose == 1
		fprintf('The oil price has changed %.2f in %s for %d gens\n', delta, year, length(find(idx)))
    end  


    % For ng
	idx = strcmp(mpc.genfuel, 'ng');
	if isfield(mpc, 'newgen')
		idx = idx & mpc.newgen ~= 1; % for non-new generator
	end
	delta = yearInfo{'ngPrice', idxYear} - yearInfo{'ngPrice', idxYear - 1};
    mpc.gencost(idx, 5) = max(mpc.gencost(idx, 5) + ...
        delta * mpc.gen_aux_data(idx, 8), mpc.gen_aux_data(idx, 9)); % 8th column is the heat rate
	if verbose == 1
		fprintf('The ng price has changed %.2f in %s for %d gens\n', delta, year, length(find(idx)))
    end

	% For ngcc, ngt and ngccccs
	idx = strcmp(mpc.genfuel, 'ngccccs') | strcmp(mpc.genfuel, 'ngcc') | strcmp(mpc.genfuel, 'ngt');
	% 8th column is the heat rate; 9th column is the VOM
	mpc.gencost(idx, 5) = yearInfo{'ngPrice', idxYear} * mpc.gen_aux_data(idx, 8) + ...
		mpc.gen_aux_data(idx, 9);

	% Add fuel price in mpc
	mpc.coal_price = yearInfo{'coalPrice', idxYear};
	mpc.ng_price = yearInfo{'ngPrice', idxYear};
	mpc.oil_price = yearInfo{'oilPrice', idxYear};

	if verbose == 1
		fprintf('The new ng price has set in %s for %d gens\n', year, length(find(idx)))
    end

	