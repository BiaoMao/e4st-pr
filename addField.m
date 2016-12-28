function mpc = addField(mpc, verbose)
%% addField: add some fields to MPC
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

	% Set default argin
	if nargin < 2
		verbose = 1; % show a little debug information
	end

    % Add number of buses and generators
    mpc.ng = size(mpc.gen, 1);
    mpc.nb = size(mpc.bus, 1);

    % Make gen built vector
    mpc.newgen = zeros(mpc.ng, 1); 

    % Debug information
    if verbose == 1
        fprintf('ng, nb, and newgen are added to mpc\n');
    end
end