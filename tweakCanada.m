function [mpc, offer] = tweakCanada(mpc, offer, caseInfo, verbose)
%% tweakCanada: make additional tweaks for Canada
%   Require capacity additions in Alberta (AB): 600 MW of wind, and 2100 MW of NGCCs;
%   Have no buildable generators in Canada, except in Alberta;
%   For generators in Canada that shut down, remove them from the dataset;
%   For generators in the US that shut down, leave them in the dataset;
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

    %% Remove buildable gen in Cananda except Alberta
    genTable = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
    genTable = join(genTable, caseInfo.locationInfo);
    idxCA = strcmp(genTable{:, 'Nation'}, 'CA') & ~strcmp(genTable{:, 'State'}, 'alberta') & ...
            ~strcmp(mpc.genfuel, 'dl') & mpc.newgen == 1;
    [mpc, offer] = removeGen(mpc, offer, idxCA);
    if verbose == 1
        fprintf('Remove %d buildable gen from Cananda except Alberta\n', length(find(idxCA)));
    end



	
	