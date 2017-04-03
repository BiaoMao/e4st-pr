function [mpc] = applyPartialRPS(mpc, caseInfo, yearInfo, year, verbose)
%% applyRPS: Apply RPS standards to partial of the system
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

        % Get Location information
        genBus = table(mpc.gen(:, 1), 'VariableNames', {'bus'});
        genBus = join(genBus, caseInfo.locationInfo);

        % Get bus idx of RE group
        if strcmp(caseInfo.rpsRE, 'all')
            idxRe = 1 : size(mpc.ng, 1);
        elseif strcmp(caseInfo.rpsRE, 'NPE')
            preExist = min(mpc.newgen);
            idxRe = mpc.newgen ~= preExist;
        elseif strcmp(caseInfo.rpsRE, 'NE-NPE')
            % Find non-preexisting gen
            preExist = min(mpc.newgen);
            idxRe = mpc.newgen ~= preExist;
            idxState = zeros(size(idxRe, 1), 1);
            % Find each NE states
            for curState = caseInfo.neRPSArea'
                idxState = idxState | strcmp(genBus{:, 'State'}, curState);
            end
            idxRe = idxRe & idxState;
        elseif strcmp(caseInfo.rpsRE, 'WECC-NPE')
            % Find non-preexisting gen
            preExist = min(mpc.newgen);
            idxRe = mpc.newgen ~= preExist;
        end

        % Get idx of load group
        if strcmp(caseInfo.rpsLoad, 'all')
            idxLoad = strcmp(mpc.genfuel, 'dl');
        elseif strcmp(caseInfo.rpsLoad, 'US')            
            idxLoad = strcmp(mpc.genfuel, 'dl') & strcmp(genBus{:, 'Nation'}, 'US');
        elseif strcmp(caseInfo.rpsLoad, 'NE')
            idxLoad = strcmp(mpc.genfuel, 'dl');
            idxState = zeros(size(idxLoad, 1), 1);
            % Find each NE states
            for curState = caseInfo.neRPSArea'
                idxState = idxState | strcmp(genBus{:, 'State'}, curState);
            end
            idxLoad = idxLoad & idxState;
        elseif strcmp(caseInfo.rpsLoad, 'WECC')
            idxLoad = strcmp(mpc.genfuel, 'dl');
        end
        
        mapRE = idxRe & (strcmp(mpc.genfuel, 'solar')...
                | strcmp(mpc.genfuel, 'wind') | strcmp(mpc.genfuel, 'oswind'));
        mapLoad = idxLoad; % total loads
        mapRPS = mapRE | mapLoad;
        % There is negative sign before dl
        coeffs = -rps * mapLoad - mapRE; % total_output <= K

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