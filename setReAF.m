function [mpc] = setReAF(mpc, fuelAf, fuelType, verbose)
%% setReAf: set availability factors for existing wind and solar
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

    % Initialize 
    nGen = size(mpc.gen, 1);
    nHours = size(fuelAf, 2);
    idx = (1 : nGen)';
	if isfield(mpc, 'newgen')   
		idx = idx & mpc.newgen ~= 1;
	end

    % Update AF by each bus 
    genIdx = find(idx & strcmp(mpc.genfuel, fuelType));
    for i = genIdx'   
    	afIdx = fuelAf{:, 'bus'} == mpc.gen(i, 1);
    	if any(afIdx)
            mpc.availability_factor(i, :) = fuelAf{afIdx, 2 : end};    		
        else
            fprintf('%s at bus %d is set with default AF. \n', fuelType, mpc.gen(i, 1));    	
        end    	
    end	

    % Debug information
    if verbose == 1
    	fprintf('Availability factors for existing %s are set.\n', fuelType);
    end