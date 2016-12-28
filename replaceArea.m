function mpc = replaceArea(mpc, caseInfo, verbose)
%% replaceArea: replace the area column in MPC
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
    busArea = caseInfo.busArea;  
    nArea = length(busArea);
    % Set the area zone     
    if(nArea == 1 && busArea == 0)
        mpc.bus(:, BUS_AREA) = 1; % for only one area
    else        
        for i = 1 : nArea - 1
            idx = mpc.bus(:, BUS_AREA) >= busArea(i) & mpc.bus(:,BUS_AREA) < busArea(i + 1);  
            mpc.bus(idx, BUS_AREA) = i;
        end
        idx = mpc.bus(:,BUS_AREA) >= busArea(i + 1);  
        mpc.bus(idx, BUS_AREA) = i + 1;
    end

    % Debug information
    if verbose == 1
        fprintf('There are %d area in mpc\n', nArea);
    end
end