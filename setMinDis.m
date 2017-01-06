function [offer] = setMinDis(mpc, offer, fuel, limit, verbose)
%% setMinDis: set minimum dispatch for certain type 
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

    iGen = strcmp(mpc.genfuel, fuel);
    offer(iGen, 13) = limit;       

    % Debug information
    if verbose == 1
        fprintf('Min dispatch of %s is set to %f\n', fuel, limit);
    end

 