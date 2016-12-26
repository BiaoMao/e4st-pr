function [mpc] = setGenAF(mpc, genAf, verbose)
%% setGenAF: set availability factors for existing generators
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

    % Initialize 
    nGen = size(mpc.gen, 1);
    nHours = size(genAf, 2);
    idx = (1 : nGen)';
	if isfield(mpc, 'newgen')   
		idx = idx & mpc.newgen ~= 1;
	end
	if ~isfield(mpc, 'availability_factor')       
    	mpc.availability_factor = ones(nGen,nHours);
    end

    % Update AF by fuel type and hour
    fuelTypes = unique(mpc.genfuel);
    for fuel = fuelTypes'
    	if strcmp(fuel, 'dl')
    		continue;
    	end
    	if ~any(strcmp(genAf.Properties.RowNames, fuel))
    		fprintf('Can not find AFs for %s.\n', fuel);
    	end
    	curIdx = idx & strcmp(mpc.genfuel, fuel);
    	mpc.availability_factor(curIdx, :) = repmat(genAf{fuel, :}, length(find(curIdx)), 1);
    end	

    % Debug information
    if verbose == 1
    	fprintf('Availability factors for existing generators are set\n')
    end