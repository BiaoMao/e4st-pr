function [mpc] = setGencost(mpc, caseInfo, verbose)
%% setGencost: Set gencost for certain fuel types
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

    % Set hydro gencost
    if isfield(caseInfo, 'gencostHydro')
        type = 'hydro';
        idx = strcmp(mpc.genfuel, type);
        mpc.gencost(idx, 5) = caseInfo.gencostHydro;
        if verbose == 1
            fprintf('Set %s gencost to %f\n', type, caseInfo.gencostHydro);
        end 
    end

    % Set refuse gencost
    if isfield(caseInfo, 'gencostRefuse')
        type = 'refuse';
        idx = strcmp(mpc.genfuel, type);
        mpc.gencost(idx, 5) = caseInfo.gencostRefuse;
        if verbose == 1
            fprintf('Set %s gencost to %f\n', type, caseInfo.gencostRefuse);
        end 
    end
end