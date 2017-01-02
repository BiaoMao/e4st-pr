function [mpc] = setHydroNonEi(mpc, caseInfo, verbose)
%% setHydroNonEi: set hydro AFs and CFs for ercot and wecc
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

    define_constants;
    groups = size(caseInfo.hydroAf, 2);
    n_hydro = size(caseInfo.hydroMap, 1);            
    map = zeros(n_hydro, size(mpc.gen, 1));
    cap = zeros(n_hydro, 1);
    coeff = ones(size(mpc.gen, 1), 1);
    coeff_type = ones(n_hydro, 1);
    idx = 1;
    for i = 1:groups
        idx_members = caseInfo.hydroMap{:,i};
        n_members = length(idx_members);
        % Set hydro CF
        cap(idx: idx+n_members-1) = mpc.gen(idx_members, PMAX) * caseInfo.hydroCf(i);
        % Re-build map for each hyro
        for j = 1:n_members
            map(idx, idx_members(j)) = 1;
            idx = idx + 1;
        end
    end
    mpc.total_output.map = map;
    mpc.total_output.cap = cap;
    mpc.total_output.coeff = coeff;
    mpc.total_output.type = coeff_type;

    % Debug information
    if verbose == 1
        fprintf('Hydro AFs and CFs are set\n');
    end

 