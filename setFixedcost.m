function [offer] = setFixedcost(mpc, offer, caseInfo, verbose)
%% setFixedcost: Set fixed cost for certain fuel types
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

    % Set hydro fixed cost
    if isfield(caseInfo, 'fixedcostHydro')
        type = 'hydro';
        idx = strcmp(mpc.genfuel, type);
        offer(idx, 1) = caseInfo.fixedcostHydro;
        if verbose == 1
            fprintf('Set %s fixed cost to %f\n', type, caseInfo.fixedcostHydro);
        end 
    end

    % Set refuse fixed cost
    if isfield(caseInfo, 'fixedcostRefuse')
        type = 'refuse';
        idx = strcmp(mpc.genfuel, type);
        offer(idx, 1) = caseInfo.fixedcostRefuse;
        if verbose == 1
            fprintf('Set %s fixed cost to %f\n', type, caseInfo.fixedcostRefuse);
        end 
    end
end