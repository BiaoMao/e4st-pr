function genData = combineGenData(mpc, offer, caseInfo)
%% combineGenData: combine generator data into a single table
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    bus = mpc.gen(:, 1);
    fuelType = mpc.genfuel;
    variableCost = mpc.gencost(:, 5);
    fixedCost = offer(:, 1);
    heatRate = mpc.gen_aux_data(:, 8);
    fuelCost = mpc.gen_aux_data(:, 7);
    vomCost = mpc.gen_aux_data(:, 9);
    genData = table(bus, fuelType, variableCost, fixedCost, fuelCost, heatRate, vomCost);
    genData = join(genData, caseInfo.locationInfo(:, {'bus', 'Longtitude', 'Latitude'}));
end