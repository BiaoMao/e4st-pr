function [out, elastyTable, res] = testDemandResponse(file1, file2, contIdx)
%% testDemandResponse: calculate actual elasticity of demand response
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
    load1 = total_load(result.(curYear).cont(contIdx).bus,...
            result.(curYear).cont(contIdx).gen, loadZone, opt);
    price1 = result.(curYear).cont(contIdx).bus(:, 14) / caseInfo.probability(contIdx + 1);

    % For result file 2
    load(file2);
    curYear = ['Y' num2str(result.year)];    
    loadZone = (1 : mpc.nb)';
    opt = struct('type', 'BOTH', 'nominal', 0);
    load2 = total_load(result.(curYear).cont(contIdx).bus,...
            result.(curYear).cont(contIdx).gen, loadZone, opt);
    price2 = result.(curYear).cont(contIdx).bus(:, 14) / caseInfo.probability(contIdx + 1);
    
    % Nominal loads
    % nArea = length(unique(mpc.bus(:, 7)));
    % loadMax = total_load(mpc.bus, mpc.gen, loadZone, opt);
    % for i = 1 : nArea
    %     idxArea = mpc.bus(:, 7) == i;
    %     loadMax(idxArea) = loadMax(idxArea) * caseInfo.loadAf{i, contIdx + 1};  
    % end     

    % Report demand curve
    idxYear = find(strcmp(yearInfo.Properties.VariableNames, ['Y', num2str(curYear)]));
    idx = strcmp(mpc.genfuel, 'dl');
    mpcCont = apply_changes(contIdx, mpc, contab); 
    busRes = getBusRes(result.Y0, caseInfo);
    gencostArr = zeros(mpc.nb, 20);
    defaultLoad = zeros(mpc.nb, 1);
    defaultPrice = zeros(mpc.nb, 1);   
    % Gencost, defaultLoad and defaultPrice
    for i = 1 : mpc.nb
        idxCur = idx & mpc.gen(:, 1) == mpc.bus(i, 1);
        if ~any(idxCur)
            continue;
        end
        gencostArr(i, :) = mpcCont.gencost(idxCur, 5 : end);
        defaultPrice(i, :) = busRes.lmpHourly(i, contIdx + 1);
        defaultLoad(i, :) = -gencostArr(i, 5); % the third step loads
    end
    loadMax = total_load(mpcCont.bus, mpcCont.gen, loadZone, opt);

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
    res.loadMax = loadMax;
    res.gencostArr = gencostArr;
    res.defaultLoad = defaultLoad;
    res.defaultPrice = defaultPrice;
    out = [table(load1, load2, price1, price2,...
        loadMax, defaultLoad, defaultPrice), array2table(gencostArr)];
end