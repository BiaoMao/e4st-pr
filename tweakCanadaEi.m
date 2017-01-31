function [mpc, offer] = tweakCanadaEi(mpc, offer, caseInfo, verbose)
%% tweakCanada: make additional tweaks for Canada EI
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

    %% Get location information
    genTable = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
    genTable = join(genTable, caseInfo.locationInfo);    

    % Remove the coal-fired generators in Ontario
    idx = strcmp(genTable{:, 'State'}, 'ontario') & strcmp(mpc.genfuel, 'coal');
    [mpc, offer] = removeGen(mpc, offer, idx);
    genTable(idx, :) = []; % consistent with dim
    if verbose == 1
        fprintf('Remove %d coal-fired generators in Ontario\n', length(find(idx)));
    end
    
    % Halve the generation capacity of the coal-fired generators in NB and NS.
    idx = (strcmp(genTable{:, 'State'}, 'new brunswick') | strcmp(genTable{:, 'State'}, 'nova scotia'))...
        & strcmp(mpc.genfuel, 'coal');
    mpc.gen(idx, 9) = mpc.gen(idx, 9) / 2;
    offer(idx, 2) = offer(idx, 2) / 2;
    if verbose == 1
        fprintf('Halve capacity of %d coal-fired generators in NB and NS\n', length(find(idx)));
    end



	
	