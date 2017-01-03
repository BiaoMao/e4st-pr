function busRes = getGenRes(mpc, offer, result, caseInfo, year)
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

    % Hourly results
    busRes.genTableHourly = zeros(size(result.base.gen, 1), caseInfo.hours);
    for i = 1 : caseInfo.hours
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
        busRes.genTableHourly(:,i) = gen(:,2) * caseInfo.hours(i);
    end

    %% Get summary results
    busRes.genTable = table(result.base.gen(:,1), result.fuel, result.newgen,...
                        'VariableNames', {'bus', 'fuel', 'newgen'});
    busArea = array2table(bus(:, [1 7]), 'VariableNames', {'bus', 'busArea'}); % bus area 
    busRes.genTable = join(busRes.genTable, busArea);
    
    % Compute generator weighted LMP
    annualGen = sum(busRes.genTableHourly, 2);
    LMPBygen = sum(genPrice .* busRes.genTableHourly, 2) ./ annualGen;
    % Re-compute zero generation without weighting
    idx = annualGen == 0;
    LMPBygen(idx, :) = genPrice(idx, :) * caseInfo.probability;

    % CO2, NOX, SO2
    totalHours = sum(caseInfo.hours);
    CO2 = result.auxdata(:,1) .* annualGen;
    NOx = result.auxdata(:,2) .* annualGen;
    SO2 = result.auxdata(:,3) .* annualGen;         
    damNOx = annualGen .* result.auxdata(:,5);
    damSO2 = annualGen .* result.auxdata(:,6);

    % Calculate CO2 damage using US and CA price
    % damCO2 = CO2 * policyInfo.damageValue(year, 1); % first column is co2 damage   
    
    % Calculate capacity
    usedCap = result.reserve.qty.Rp_pos;
    isNewgen = result.newgen == 1;
    invests = usedCap;
    invests(~isNewgen) = 0;

    % Calculate retirements
    % Exclude Year 0 and some fuel types in the policyInfo
    retires = result.cap - usedCap;
    retires(isNewgen) = 0;
    if year == 1 || year == 0
        retires(:) = 0;
    end
    
    fixedCost = result.fixedcost .* result.reserve.qty.Rp_pos * totalHours;
    variableCost = result.gencost .* annualGen;

    % Calculate tax and insurance: for used cap only
    fuels= busRes.genTable{:,2};
    fuelTypes = unique(fuels);
    tax = zeros(size(result.base.gen,1), 1);
    insurance = zeros(size(result.base.gen,1), 1);
    for i = 1: length(fuelTypes)
        idxGen = find(strcmp(fuels,fuelTypes{i}));
        idxInfo = find(strcmp(caseInfo.fuelName,fuelTypes{i}));
        tax(idxGen) = result.reserve.qty.Rp_pos(idxGen) * caseInfo.genInfo(idxInfo,5) *totalHours;
        insurance(idxGen) =result.reserve.qty.Rp_pos(idxGen) * caseInfo.genInfo(idxInfo,6) *totalHours;
    end

    % Combine table
    busRes.genTable = [busRes.genTable table(annualGen, usedCap, retires, invests, fixedCost, variableCost,...
                    tax, insurance, CO2, NOx, SO2, damCO2, damNOx, damSO2, LMPBygen)];   

    % Set the dl values to zero
    idxDl = strcmp(busRes.genTable{:,'fuel'}, 'dl');
    busRes.genTable{idxDl, 5:end} = 0;      

    %% Add addtional and optional part
    if isfield(caseInfo, 'busLoc')
        busRes.genTable = join(busRes.genTable, caseInfo.busLoc);
    end 
end