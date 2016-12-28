function mpc = relaxMPC(mpc, verbose)
%% relaxMPC: relax some constraints in MPC 
%   Set RAMP_10 to Inf (for all generators) to ensure that the physical ramp rate never 
%   limits the change in dispatch of a unit from one hour type to another;
%   Set PMIN to zero;
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

    define_constants;
    mpc.gen(:, RAMP_10) = Inf;
    iPmin = mpc.gen(:, PMIN) > 0;
    iGen = ~strcmp(mpc.genfuel, 'dl');
    mpc.gen(iPmin & iGen, PMIN) = 0;   

    % Debug information
    if verbose == 1
        fprintf('RAMP_10 is set to Inf; PMIN is set to zero\n');
    end
end