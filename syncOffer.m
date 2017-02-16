function offer = syncOffer(mpc, offer, verbose)
%% syncOffer: make offer(:, 2) = mpc.gen(:, PMAX)
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
    
    offer(:, 2) = max(mpc.gen(:, PMAX), -mpc.gen(:, PMIN));

    % Debug information
    if verbose == 1
        fprintf('Sync offer with PMAX\n');
    end
end