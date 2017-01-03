function mpc = addCO2Cap(mpc, isInCap, cap, verbose)
%% addCO2Cap: add CO2 cap to MPC
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

	% Set default argin
	if nargin < 5
		verbose = 1; % show a little debug information
	end

    % Check if total_output filed is initialized
    if ~isfield(mpc, {'total_output'})
        mpc.total_output.map = [];
        mpc.total_output.cap = [];
        mpc.total_output.coeff = [];
        mpc.total_output.type = [];
    end

    % Get bus idx of non rggi
    busGIS = array2table([mpc.bus(:, 1), isInCap], 'VariableNames', {'bus', 'map'});
    genBus = table(mpc.gen(:, 1), 'VariableNames',{'bus'});
    mapTable = join(genBus, busGIS);
    cap_map = mapTable{:, 'map'};

    mpc.total_output.map = [mpc.total_output.map; cap_map'];
    mpc.total_output.cap = [mpc.total_output.cap; cap]; 

    % Second column of coeff is CO2 emission rate
    if size(mpc.total_output.coeff, 2) == 1
        mpc.total_output.coeff(:, 2) = mpc.gen_aux_data(:,1);
    end

    % select the second columns of coeffs, first row is ones
    mpc.total_output.type = [mpc.total_output.type; 2]; 

    % Debug information
    if verbose == 1
        fprintf('CO2 cap is added\n');
    end
end