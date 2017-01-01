function [mpc] = setReAF(mpc, fuelInfo, fuelType, idxGen, verbose)
%% setReAf: update availability factors for wind and solar
%	AFs should be set by default values first	
%	idxGen is optional: only updates parts of the AFs
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

	% Set default argin
	if nargin < 5
		verbose = 1; % show a little debug information
	end

    % Initialize 
    nGen = size(mpc.gen, 1);
    startCol = find(strcmp(fuelInfo.Properties.VariableNames, 'C1')); % Starting of AF table
    idx = (1 : nGen)';
    if exist('idxGen', 'var')   
		idx(~idxGen) = 0;
	end

    % Update AF by each bus 
    genIdx = find(idx & strcmp(mpc.genfuel, fuelType));
    for i = genIdx'   
    	afIdx = fuelInfo{:, 'bus'} == mpc.gen(i, 1);
    	if any(afIdx)
            mpc.availability_factor(i, :) = fuelInfo{afIdx, startCol : end};    		
        else
            fprintf('%s at bus %d is set with default AF\n', fuelType, mpc.gen(i, 1));    	
        end    	
    end	

    % Debug information
    if verbose == 1
    	fprintf('Availability factors for existing %s are set\n', fuelType);
    end