function [mpc] = addCO2Tax(mpc, yearInfo, year, verbose)
%% addCO2Tax: add CO2 tax
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

    curYear = ['Y' num2str(year)];
    idx = find(strcmp(yearInfo.Properties.VariableNames, curYear));

    mpc.gencost(:, 5) = mpc.gencost(:, 5) + mpc.gen_aux_data(:, 1) * yearInfo{'emissPriceCO2', idx};

	if verbose == 1
    	fprintf('Curren CO2 tax is %f\n', yearInfo{'emissPriceCO2', idx});
    end  
	