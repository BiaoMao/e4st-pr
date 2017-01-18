function [mpc, contab] = makeNewDr(mpc, contab, caseInfo, yearInfo, busRes, year, mode, verbose)
%% makeNewDr: make the step gencost of the demand response for base hour and contingency hours
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    % Set default argin
    if nargin < 8
        verbose = 1; % show a little debug information
    end

    % Extend mpc columns if necessary
	if size(mpc.gencost, 2) == 6
        % 10 steps of pairs; first pair is already included
        mpc.gencost = [mpc.gencost zeros(mpc.ng, 18)];
    end

    % Make step gencost for base hour
    idxDl = strcmp(mpc.genfuel, 'dl');    
    mpc.gencost(idxDl, :) = makeDrStep(mpc, caseInfo, yearInfo, busRes, 1, year, mode, verbose);

    % Make step gencost of contab for contingency hours
    contabDr = makeContDr(mpc, caseInfo, yearInfo, busRes, year, mode, verbose);
    contab = [contab; contabDr];

    % Debug information
    if verbose == 1
        fprintf('Demand responses in Year %d are set\n', year);
    end
end