function [elastyTable, res] = testDRBase(file1, hour)
%% testDRBase: calculate actual elasticity of demand response
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    % For result file 1
    load(file1);
    curYear = ['Y' num2str(result.year)];    
    loadZone = (1 : mpc.nb)';
    opt = struct('type', 'BOTH', 'nominal', 0);
    load1 = total_load(result.(curYear).cont(hour).bus,...
            result.(curYear).cont(hour).gen, loadZone, opt);
    price1 = result.(curYear).cont(hour).bus(:, 14) / caseInfo.probability(hour);

    % For result file 2
    idxYear = find(strcmp(yearInfo.Properties.VariableNames, curYear));
    baseYear = 'Y0';    
    loadZone = (1 : mpc.nb)';
    opt = struct('type', 'BOTH', 'nominal', 0);
    load2 = total_load(result.(baseYear).cont(hour).bus,...
            result.(baseYear).cont(hour).gen, loadZone, opt);
    nArea = length(unique(mpc.bus(:, 7)));
    for i = 1 : nArea
        idxArea = loadArea == i;
        load2(idxArea) = load2(idxArea) * caseInfo.loadGrowth(i).^sum(yearInfo{'loadYearDelta', 2 : idxYear});  
    end 
    price2 = result.(baseYear).cont(hour).bus(:, 14) / caseInfo.probability(hour);

    % Calculate elasticity
    elastyTable = log(load1 ./ load2) ./ log((price1 + 60) ./ (price2 + 60));
    res.load1 = load1;
    res.load2 = load2;
    res.price1 = price1;
    res.price2 = price2;
    
    res.totalLoad1 = sum(load1);
    res.totalLoad2 = sum(load2);
    res.lmp1 = price1' * load1 / res.totalLoad1;
    res.lmp2 = price2' * load2 / res.totalLoad2;
    res.totalElasty = log(res.totalLoad1 / res.totalLoad2) /...
                      log((res.lmp1 + 60) / (res.lmp2 + 60));
end