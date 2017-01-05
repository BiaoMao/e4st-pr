function [mpc] = addInterface(mpc, caseInfo, verbose)
%% addInterface: add interface constraints
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

	% Set default argin
	if nargin < 4
		verbose = 1; % show a little debug information
	end

    mpc.if.map = caseInfo.interfaceMap{:, :};
    % Add 4850 MW power flow interface in NY
    mpc.if.lims = caseInfo.interfaceLims{:, :};
    mpc = toggle_iflims(mpc, 'on');

	if verbose == 1
    	fprintf('Interface constraint is added\n')
    end  
	