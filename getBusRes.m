function busRes = getBusRes(result, caseInfo)
%% getBusRes: summarrize the results by bus
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

	%% Get hourly results
    loadHourly = zeros(size(result.base.bus, 1), caseInfo.nHours);
    lmpHourly = zeros(size(result.base.bus, 1), caseInfo.nHours);
    for i = 1: caseInfo.nHours
        if i == 1
            bus = result.base.bus;
            gen = result.base.gen;
        else
            bus = result.cont(i - 1).bus;
            gen = result.cont(i - 1).gen;
        end     
        loadHourly(:,i) = total_load(bus, gen, (1: size(bus, 1))') * caseInfo.hours(i);
        lmpHourly(:,i) = bus(:,14) / caseInfo.probability(i);
    end
    % Order of index is exactly the same with mpc.bus
    busRes.loadHourly = loadHourly;
    busRes.lmpHourly = lmpHourly;

    %% Get summary results
    annualLoads = sum(loadHourly, 2);
    annualPrices = sum(loadHourly .* lmpHourly, 2) ./ annualLoads;
    
    % Re-compute zero loads lmp without load weighting
    idx = annualLoads == 0;
    annualPrices(idx) = lmpHourly(idx,:) * caseInfo.probability;

    %% Combine table
    busArea = bus(:, 7); % bus area 
    bus = bus(:, 1); % in order to set            
    busRes.busTable = table(bus, busArea, annualLoads, annualPrices);

    %% For optional and additional parts
    if nargin >= 3
        busRes.busTable = join(busRes.busTable, locInfo);
    end

    %% Calclate system result
    busRes.sysLoad = sum(annualLoads);
    busRes.sysLMP = annualPrices' * annualLoads / busRes.sysLoad;
end