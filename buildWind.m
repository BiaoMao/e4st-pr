function [mpc, offer] = buildWind(mpc, offer, caseInfo, verbose)
%% buildWind: add buildable wind
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

    % Initial data
    if isfield(caseInfo, 'newWindInfo')
        % If there are two sets of windInfo
        windInfo = caseInfo.newWindInfo;
    else
        windInfo = caseInfo.windInfo;
    end
    genInfo = caseInfo.genInfo;   

    % Build wind
    fuelType = 'wind';
    iGeninfo = strcmp(genInfo.Properties.RowNames, fuelType);
    busId = windInfo{:, 'bus'}; % get unique bus id
    numBus = length(busId);

    % Initial all the new table
    newGen = zeros(numBus,21);
    newGencost = zeros(numBus, size(mpc.gencost, 2));
    newGenfuel = cell(numBus,1);
    newGenaux = zeros(numBus,10);
    newGenindex = zeros(numBus,1);
    newOffer = zeros(numBus,13);

    % Add new gen table             
    newGen(:, 1) = busId; % bus number
    newGen(:,6) = 1; % voltage magnitude setpoint
    newGen(:,7) = 100; % MVA base
    newGen(:,8) = 1; % gen status
    newGen(:,9) = Inf; % PMAX
    newGen(:,18) = Inf; % ramp rate for 10 min 

    % Add new gencost table
    newGencost(:, 1:4) = repmat([2 0 0 2], numBus, 1);
    newGencost(:, 5) = genInfo{iGeninfo, 'Gencost'}; % gencost

    % Add new gen fuel table
    newGenfuel(:, 1) = {fuelType};

    % Add new gen aux data table
    newGenaux = repmat(genInfo{iGeninfo, 8 : 17}, numBus, 1);

    % Add new genbuilt and newgen
    newGenindex = ones(numBus, 1);

    % Add new AFs
    startCol = find(strcmp(windInfo.Properties.VariableNames, 'C1')); % Starting of AF table
    newAf = windInfo{:, startCol : end};

    % Add new offer table
    newOffer(:, 1) = sum(genInfo{iGeninfo, {'Cost2Keep', 'Cost2Build', 'Tax', 'Insurance'}}, 2); % fixed cost
    newOffer(:, 2) = windInfo{:, 'Cap'}; % Installed Cap
    newOffer(:, 3) = 0;
    newOffer(:, 4) = Inf;

    % Update in the mpc and offer            
    mpc.gen = [mpc.gen; newGen];
    mpc.gencost = [mpc.gencost; newGencost];
    mpc.genfuel = [mpc.genfuel; newGenfuel];
    mpc.gen_aux_data = [mpc.gen_aux_data; newGenaux];          
    mpc.newgen = [mpc.newgen; newGenindex];
    mpc.availability_factor = [mpc.availability_factor; newAf];
    mpc.ng = size(mpc.gen,1);
    mpc.nb = size(mpc.bus,1);
    offer = [offer; newOffer];

    % Debug information
    if verbose == 1
        fprintf('Buildable %s are added\n', fuelType);
    end
end