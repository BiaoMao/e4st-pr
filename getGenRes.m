function genRes = getGenRes(mpc, offer, result, caseInfo, yearInfo, year)
%% getGenRes: summarrize the results by gen
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    % Initial input data
    locInfo = caseInfo.locationInfo;
	year = ['Y', num2str(year)];
	idxYear = find(strcmp(yearInfo.Properties.VariableNames, year));

    % Hourly results
    genRes.genTableHourly = zeros(size(result.base.gen, 1), caseInfo.nHours);
    genPrice = zeros(size(result.base.gen, 1), caseInfo.nHours);
    for i = 1 : caseInfo.nHours
        if i == 1
            gen = result.base.gen;
            bus = result.base.bus;                    
        else
            gen = result.cont(i - 1).gen;
            bus = result.cont(i - 1).bus;
        end     
        genBus = table(gen(:,1), 'VariableNames',{'bus'});
        busPrice = table(bus(:,1), bus(:,14), 'VariableNames',{'bus','price'});
        genLmp = join(genBus, busPrice);
        genPrice(:, i) = genLmp{:,2} / caseInfo.probability(i); % LMP at each generator bus
        genRes.genTableHourly(:,i) = gen(:,2) * caseInfo.hours(i);
    end
    genRes.lmpHourly = genPrice;

    %% Get summary results
    genRes.genTable = table(result.base.gen(:,1), mpc.genfuel, mpc.newgen,...
                        'VariableNames', {'bus', 'fuel', 'newgen'});
    busArea = array2table(bus(:, [1 7]), 'VariableNames', {'bus', 'busArea'}); % bus area 
    genRes.genTable = join(genRes.genTable, busArea);
    
    % Compute generator weighted LMP
    annualGen = sum(genRes.genTableHourly, 2);
    LMPBygen = sum(genPrice .* genRes.genTableHourly, 2) ./ annualGen;
    % Re-compute zero generation without weighting
    idx = annualGen == 0;
    LMPBygen(idx, :) = genPrice(idx, :) * caseInfo.probability;

    % CO2, NOX, SO2
    caseInfo.nHours = sum(caseInfo.hours);
    CO2 = mpc.gen_aux_data(:,1) .* annualGen;
    NOx = mpc.gen_aux_data(:,2) .* annualGen;
    SO2 = mpc.gen_aux_data(:,3) .* annualGen;    

    % Calculate damages     
    isNewgen = mpc.newgen == 1;
    damNOx(~isNewgen, 1) = annualGen(~isNewgen, 1) .* mpc.gen_aux_data(~isNewgen, 5) * yearInfo{'damageNOxPE', idxYear};
    damSO2(~isNewgen, 1) = annualGen(~isNewgen, 1) .* mpc.gen_aux_data(~isNewgen, 6) * yearInfo{'damageSO2PE', idxYear};
    damNOx(isNewgen, 1) = annualGen(isNewgen, 1) .* mpc.gen_aux_data(isNewgen, 5) * yearInfo{'damageNOxNPE', idxYear};
    damSO2(isNewgen, 1) = annualGen(isNewgen, 1) .* mpc.gen_aux_data(isNewgen, 6) * yearInfo{'damageSO2NPE', idxYear};

    % Calculate CO2 damage using US and CA price
    damCO2 = CO2 * yearInfo{'damageCO2', idxYear}; % first column is co2 damage   
    
    % Calculate capacity
    usedCap = result.reserve.qty.Rp_pos;
    assignedCap = offer(:, 2);
    investCap = usedCap;
    investCap(~isNewgen) = 0;

    % Calculate retirements
    % Exclude Year 0 and some fuel types in the policyInfo
    shutDownCap = assignedCap - usedCap;
    shutDownCap(isNewgen) = 0;
    % if idxYear == 1
    %     shutDownCap(:) = 0;
    % end
    
    % Calculate costs    
    fixedCost = offer(:, 1) .* usedCap * caseInfo.nHours;
    variableCost = mpc.gencost(:, 5) .* annualGen;

    % Calculate tax and insurance: for used cap only
    fuels = genRes.genTable{:, 2};
    fuelTypes = unique(fuels);
    tax = zeros(size(result.base.gen,1), 1);
    insurance = zeros(size(result.base.gen,1), 1);
    for i = 1: length(fuelTypes)
        idxGen = find(strcmp(fuels, fuelTypes{i}));
        if strcmp(fuelTypes{i}, 'oswind')
            tax(idxGen) = 0;
            insurance(idxGen) = 0;
            continue;
        end
        tax(idxGen) = usedCap(idxGen) * caseInfo.genInfo{fuelTypes{i}, 'Tax'} * caseInfo.nHours;
        insurance(idxGen) = usedCap(idxGen) * caseInfo.genInfo{fuelTypes{i}, 'Insurance'} * caseInfo.nHours;
    end

    % Combine table
    genRes.genTable = [genRes.genTable table(annualGen, usedCap,    shutDownCap, investCap, fixedCost, variableCost,...
                    tax, insurance, CO2, NOx, SO2, damCO2, damNOx, damSO2, LMPBygen)];   

    % Set the dl values to zero
    idxDl = strcmp(genRes.genTable{:,'fuel'}, 'dl');
    genRes.genTable{idxDl, 5 : end} = 0;      

    %% Add addtional and optional part
    if isfield(caseInfo, 'locationInfo')
        genRes.genTable = join(genRes.genTable, caseInfo.locationInfo);
    end 
end