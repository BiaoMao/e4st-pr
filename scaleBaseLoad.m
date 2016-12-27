function [mpc] = scaleBaseLoad(mpc, scalings, verbose)
%% scaleAreaLoad: scale loads in a area with scalling in base hour
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
    BUS_AREA = 7;        
    nArea = length(scalings);
    if nArea == 1 % there is only one area
        [mpc.bus,mpc.gen] = scale_load(scalings, mpc.bus, mpc.gen);
    else % there are multi-area
        loadZone = mpc.bus(:, BUS_AREA);
        opt = struct('pq','P');
        [mpc.bus, mpc.gen] = scale_load(scalings, mpc.bus, mpc.gen, loadZone, opt);
    end 

    % Debug information
    if verbose == 1
        fprintf('Base loads are scaled by ');
        fprintf('%.2f\t', scalings);
        fprintf('\n');
    end

 