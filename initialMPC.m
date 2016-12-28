function mpc = initialMPC(mpc, caseInfo, verbose)
%% initialMPC: preprocess MPC
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

	mpc = selIsland(mpc, caseInfo, verbose);
	mpc = relaxMPC(mpc, verbose);
	mpc = replaceArea(mpc, caseInfo, verbose);
	mpc = addField(mpc, verbose);
end