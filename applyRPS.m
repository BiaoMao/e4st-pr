function [mpc] = applyRPS(mpc, yearInfo, year, verbose)
%% applyRPS: Apply RPS standards
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

    % Format the year as a string
	year = ['Y', num2str(year)];
	idxYear = find(strcmp(yearInfo.Properties.VariableNames, year));  

    % Check if total_output filed is initialized
    if ~isfield(mpc, 'total_output') || isempty(mpc.total_output)
        mpc.total_output.map = [];
        mpc.total_output.cap = [];
        mpc.total_output.coeff = [];
        mpc.total_output.type = [];
    end

    % Apply RPS to non pre-existing solar and wind generators
    if any(strcmp(yearInfo.Properties.RowNames, 'reRPS'))   
        rps = yearInfo{'reRPS', idxYear};

        % Get bus idx of non pre-existing 
        preExist = min(mpc.newgen);
        mapRE = (mpc.newgen ~= preExist) & (strcmp(mpc.genfuel, 'solar')...
                | strcmp(mpc.genfuel, 'wind') | strcmp(mpc.genfuel, 'oswind'));
        mapLoad = ~(strcmp(mpc.genfuel, 'dl')); % total generations
        mapRPS = mapRE | mapLoad;
        coeffs = rps * mapLoad - mapRE; % total_output <= K

        mpc.total_output.map = [mpc.total_output.map; mapRPS'];
        mpc.total_output.cap = [mpc.total_output.cap; 0]; 
        mpc.total_output.coeff = [mpc.total_output.coeff coeffs];
        % select the latest columns of coeffs, first row is ones
        mpc.total_output.type = [mpc.total_output.type; size(mpc.total_output.coeff, 2)]; 

        if verbose == 1
            fprintf('Apply %f RPS to non pre-existing wind and solar\n', rps);
        end 
    end
end