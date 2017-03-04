function [mpc, offer, caseInfo] = buildOSWind(mpc, offer, caseInfo, location, oswSize, verbose)
%% buildOSWind: add buildable offshore wind
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    % Set default argin
    if nargin < 6
        verbose = 1; % show a little debug information
    end

    % Initial data
    windInfo = caseInfo.osWind;
    idxNew = strcmp(windInfo{:, 'location'}, location) &...
            strcmp(windInfo{:, 'size'}, oswSize);
    fprintf('Curent osw is %s and %s\n', location, oswSize);
    newCap = windInfo{idxNew, 'cap'}(1);       
    preIdx = size(mpc.gen, 1);

    % Build wind
    fuelType = 'oswind';
    busId = windInfo{idxNew, 'bus'}; % get unique bus id
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
    newGen(:,9) = newCap; % PMAX
    newGen(:,18) = Inf; % ramp rate for 10 min 

    % Add new gencost table
    newGencost(:, 1:4) = repmat([2 0 0 2], numBus, 1);
    newGencost(:, 5) = windInfo{idxNew, 'gencost'}; % gencost

    % Add new gen fuel table
    newGenfuel(:, 1) = {fuelType};

    % Add new gen aux data table
    % All zeros for wind

    % Add new genbuilt and newgen
    newGenindex = ones(numBus, 1);

    % Add new AFs
    startCol = find(strcmp(windInfo.Properties.VariableNames, 'C1')); % Starting of AF table
    newAf = windInfo{idxNew, startCol : end};

    % Add new offer table
    % Fixed cost, cost2k, tax and insurance are zero        
    newOffer(:, 1) = windInfo{idxNew, 'cost2build'}; 
    newOffer(:, 2) = newCap; % Installed Cap
    newOffer(:, 3) = 0;
    newOffer(:, 4) = Inf;
    % Update the genInfo table
    caseInfo.genInfo{fuelType, 'Cost2Build'} = newOffer(1, 1);
    caseInfo.genInfo{fuelType, 'InstallCap'} = newOffer(1, 2) ;

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

    % Add capacity constraints
    if ~isfield(mpc, 'caplim') || isempty(mpc.caplim)
        mpc.caplim.map = [];
        mpc.caplim.max = [];
        mpc.caplim.min = [];
    end
    idxOsw = preIdx + 1 : preIdx + numBus;
    idxGen = zeros(1, size(mpc.gen, 1));
    idxGen(idxOsw) = 1;
    mpc.caplim.map = [mpc.caplim.map; idxGen];
    mpc.caplim.max = [mpc.caplim.max; newCap];
    mpc.caplim.min = [mpc.caplim.min; newCap];

    % Debug information
    if verbose == 1
        fprintf('Buildable %s %s are added at %s\n', oswSize, fuelType, location);
    end
end