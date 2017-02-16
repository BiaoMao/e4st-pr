classdef E4STResult < matlab.mixin.Copyable
%E4ST results class to summarize the results
%   Create table results for E4ST

%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
   
    properties (Constant)

    end
    
        %% sumBynation: summary result by US and CA
    methods (Static)
        function [outputs] = sumBynation(genTable, busTable)
            %% Get summary by nation       

            % Remove 'dl' row
            idxDl = strcmp(genTable.genTable{:, 'fuel'}, 'dl');            
            genTable.genTable(idxDl, :) = [];    

            func = @sum;
            % Only pass the variables that you want to sum
            inputVariables = {'annualGen','usedCap', 'shutDownCap', 'investCap', 'fixedCost',...
                            'variableCost', 'tax', 'insurance', 'CO2', 'NOx', 'SO2',...
                            'damCO2', 'damNOx', 'damSO2'};
            outputs = varfun(func, genTable.genTable,...
                    'GroupingVariables','nation', 'InputVariables', inputVariables);
            LMP = E4STResult.weightPrice(busTable.busTable(:, {'nation','annualLoads','LMPToLoad'}));
            genPrice = E4STResult.weightPrice(genTable.genTable(:, {'nation','annualGen','LMPToGen'}));
            outputs = [outputs table(LMP) table(genPrice)];
        end


        %% sumBystate: summary result by states
        function [outputs] = sumBystate(genTable, busTable)
            %% Get summary by state       

            % Remove 'dl' row
            idxDl = strcmp(genTable.genTable{:, 'fuel'}, 'dl');            
            genTable.genTable(idxDl, :) = [];  

            func = @sum;
            % Only pass the variables that you want to sum
            inputVariables = {'annualGen','usedCap', 'shutDownCap', 'investCap', 'fixedCost',...
                            'variableCost', 'tax', 'insurance', 'CO2', 'NOx', 'SO2',...
                            'damCO2', 'damNOx', 'damSO2'};
            outputs = varfun(func, genTable.genTable,...
                    'GroupingVariables','state', 'InputVariables', inputVariables);
            LMP = E4STResult.weightPrice(busTable.busTable(:, {'state','annualLoads','LMPToLoad'}));
            genPrice = E4STResult.weightPrice(genTable.genTable(:, {'state','annualGen','LMPToGen'}));
            outputs = [outputs table(LMP) table(genPrice)];
        end

        %% sumByfuel: summary result by fuels
        function [outputs] = sumByfuel(genRes)
            %% Get summary by fuels         
            func = @sum;
            % Only pass the variables that you want to sum
            inputVariables = {'annualGen','usedCap', 'shutDownCap', 'investCap'};
            outputs = varfun(func, genRes.genTable,...
                    'GroupingVariables', {'fuel'}, 'InputVariables', inputVariables);
            capacityFactor = outputs{:, 'sum_annualGen'} ./ (outputs{:, 'sum_usedCap'} * 8760);
            outputs = [outputs array2table(capacityFactor)];
            % Remove 'dl' row
            idxDl = strcmp(outputs{:, 'fuel'}, 'dl');
            outputs(idxDl, :) = []; 
        end

        function [outputs] = fuelRes2Row(fuelRes, fuelType, colIdx)
            %% Convert fuel table to one row        
            outputs = table;
            colSize = 5;
            for i = 1 : size(fuelRes, 1)
                if any(strcmp(fuelType, fuelRes{i, 'fuel'}))
                    % Skip fuel types 
                    continue;
                end
                curTable = [fuelRes(i, {'sum_annualGen',...
                    'sum_usedCap', 'sum_shutDownCap', 'sum_investCap', 'capacityFactor'})];
                curFuel = fuelRes{i, 'fuel'}{:};
                curTable.Properties.VariableNames = {[curFuel '_annualGen'],...
                    [curFuel '_usedCap'], [curFuel '_shutDownCap'],...
                     [curFuel '_investCap'], [curFuel '_capacityFactor']};
                curTable = curTable(:, colIdx);
                outputs = [outputs curTable];
            end
        end

        %% getMultipliers: get multipliers of from output constraints    
        %   idx: index of multipliers in cap         
        function [multipliers] = getMultipliers(result, idx)
            multipliers = result.total_output.mu(idx);  
        end

        %% sumByall: summary results for the whole system
        function [outputs] = sumBysys(genRes, busRes, caseInfo, yearInfo, result, year)
            % Only pass the variables that you want to sum

            % Remove 'dl' row
            idxDl = strcmp(genRes.genTable{:, 'fuel'}, 'dl');            
            genRes.genTable(idxDl, :) = [];  

            func = @sum;
            inputVariables = {'annualGen', 'usedCap','investCap', 'shutDownCap', ...
                            'variableCost', 'tax', 'insurance', 'costToKeep', 'objCAPEX',...
                            'CO2', 'NOx', 'SO2', 'damCO2', 'damNOx', 'damSO2'};

            %% Get gen summary for all
            genSum = varfun(func, genRes.genTable, 'InputVariables', inputVariables);
            avgLMPtoGen = genRes.genTable{:,'annualGen'}' * genRes.genTable{:,'LMPToGen'} / genSum{1, 'sum_annualGen'};
            avgLMPtoLoad = busRes.busTable{:, 'annualLoads'}' * busRes.busTable{:, 'LMPToLoad'} / genSum{1, 'sum_annualGen'};    
        
            objValue = -result.opf_results.f * sum(caseInfo.hours);

            % Calculate wind subsidy
            curYear = ['Y' num2str(year)];
            idxWind = strcmp(genRes.genTable{:, 'fuel'}, 'wind');
            totalWindSubsidy = sum(genRes.genTable{idxWind,'annualGen'}) * yearInfo{'windSub', curYear};

            % Calculate solar subsidy: 30% cost2build
            idxSolar = strcmp(genRes.genTable{:, 'fuel'}, 'solar');
            objSolarSubsidy = sum(genRes.genTable{idxSolar,'objCAPEX'}) * 0.3 / 0.7;

            % Output the summary table
            outputs = [genSum table(avgLMPtoGen, avgLMPtoLoad, totalWindSubsidy, objSolarSubsidy, objValue)];
        end

        %% sumSurplusEMFWecc: sum surplus for EMF case
        function [outputs] = sumSurplusEMFWindWecc(sysRes, multipliers, qty_Multipliers, mpc, result, caseInfo)
            % Calculate surplus
            % original obj value is hourly base

            % Map the result table to variable
            avgLMPtoLoad = sysRes{1, 'avgLMPtoLoad'};
            avgLMPtoGen = sysRes{1, 'avgLMPtoGen'};

            % Calculate CO2 payments based on cap-trade program
            numCap = length(mpc.total_output.cap);
            numMul = size(multipliers, 2);
            co2PayToGov = multipliers{1, :} * qty_Multipliers;
            co2PayFromProducer = co2PayToGov;

            govRevenue = sysRes{1,'sum_tax'} + co2PayToGov;
            % welProProfit = avgLMPtoGen * sysRes{1,'sum_annualGen'} - sysRes{1,'sum_fixedCost'} - sysRes{1,'sum_variableCost'}...
            %                 - co2PayFromProducer;
            welProProfit = avgLMPtoGen * sysRes{1,'sum_annualGen'} - sysRes{1,'sum_directCost'}...
                            - co2PayFromProducer;

            % consumerSurp = -result.opf_results.f * sum(caseInfo.nHours) + sysRes{1,'sum_fixedCost'} + sysRes{1,'sum_variableCost'}...
            %         - avgLMPtoLoad * sysRes{1,'sum_annualGen'};
            consumerSurp = -result.opf_results.f * sum(caseInfo.hours) + sysRes{1,'sum_directCost'}...
                    - avgLMPtoLoad * sysRes{1,'sum_annualGen'};
            envirSurp = -sum(sysRes{1, {'sum_damCO2', 'sum_damNOx', 'sum_damSO2'}});             
            merchanSurplus = sysRes{1,'sum_annualGen'} * (avgLMPtoLoad - avgLMPtoGen);
            totalSurplus = sum([consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus]);

            % Output the summary table
            outputs = table(consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus, totalSurplus);
        end

        %% sumSurplusBenWindWecc: sum surplus for BenWind case
        function [outputs] = sumSurplusBenWindWecc(sysRes, multipliers, mpc, result, caseInfo)
            % Calculate surplus
            % original obj value is hourly base

            % Map the result table to variable
            avgLMPtoLoad = sysRes{1, 'avgLMPtoLoad'};
            avgLMPtoGen = sysRes{1, 'avgLMPtoGen'};

            % Calculate CO2 payments based on cap-trade program
            numCap = length(mpc.total_output.cap);
            co2PayToGov = multipliers{1, 1 : 2} * result.total_output.qty(numCap - 1 : numCap) * sum(caseInfo.hours);
            co2PayFromProducer = co2PayToGov;

            % Full CAPEX amount, we need to DIVIDE by 0.14902949
            fullCAPEX = sysRes{1, 'sum_objCAPEX'} / 0.14902949;
            fullSolarSubsidy = sysRes{1, 'objSolarSubsidy'} / 0.14902949;

            % Full CAPEX amount divided by number of years
            welCAPEX = fullCAPEX / 5;
            welSolarSubsidy = fullSolarSubsidy / 5;

            producerCost = sysRes{1, 'sum_variableCost'} + welCAPEX +...
                sysRes{1, 'sum_costToKeep'} + sysRes{1, 'sum_tax'} +...
                sysRes{1, 'sum_insurance'} + co2PayFromProducer;
                   

            govRevenue = sysRes{1,'sum_tax'} + co2PayToGov -...
                        sysRes{1,'totalWindSubsidy'} - welSolarSubsidy;
            
            objProProfit = avgLMPtoGen * sysRes{1,'sum_annualGen'} -...
                (sysRes{1, 'sum_variableCost'} + sysRes{1, 'sum_objCAPEX'} +...
                sysRes{1, 'sum_costToKeep'} + sysRes{1, 'sum_tax'} +...
                sysRes{1, 'sum_insurance'});

            welProProfit = avgLMPtoGen * sysRes{1,'sum_annualGen'} - producerCost;

            envirSurp = -sum(sysRes{1, {'sum_damCO2', 'sum_damNOx', 'sum_damSO2'}});   

            merchanSurplus = sysRes{1,'sum_annualGen'} * (avgLMPtoLoad - avgLMPtoGen);

            consumerSurp = sysRes{1,'objValue'} - objProProfit - merchanSurplus;

            totalSurplus = sum([consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus]);

            proCostWOTransfer = producerCost - govRevenue;

            % Output the summary table
            outputs = table(consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus, totalSurplus,...
                producerCost, proCostWOTransfer, welCAPEX, fullCAPEX, objProProfit, co2PayToGov, welSolarSubsidy);
        end

        %% sumSurplusRPSRGGI: sum surplus for rpsrggi  case
        function [outputs] = sumSurplusRPSRGGI(genSum, genTable, busTable, caseInfo, result, yearInfo, iYear)
            % Calculate surplus
            % original obj value is hourly base

            % Map the result table to variable
            avgLMPtoLoad = genSum{1, 'avgLMPtoLoad'};
            avgLMPtoGen = genSum{1, 'avgLMPtoGen'};

            rggiCo2 = varfun(@sum, genTable.genTable, 'InputVariables', 'CO2', 'GroupingVariables', 'isRGGI');
            idxRGGI = rggiCo2{:, 'isRGGI'} == 1;

            % Calculate CO2 payments based on CO2 emission price or cost of cap-trade program
            if iYear == 1
                co2PayToGov = yearInfo.co2PriceRGGI(iYear) * rggiCo2{idxRGGI, 'sum_CO2'};
                co2PayFromProducer = 0; % already included in the gencost
            else
                multipliers = E4STResult.getMultipliers(result);
                co2PayToGov = multipliers{1, 1} * rggiCo2{~idxRGGI, 'sum_CO2'} + multipliers{1, 2} * rggiCo2{idxRGGI, 'sum_CO2'};
                co2PayFromProducer = co2PayToGov;
            end

            govRevenue = genSum{1,'sum_tax'} + co2PayToGov;
            welProProfit = avgLMPtoGen * genSum{1,'sum_annualGen'} - genSum{1,'sum_fixedCost'} - genSum{1,'sum_variableCost'}...
                            - co2PayFromProducer;

            consumerSurp = -result.f * sum(caseInfo.hours) + genSum{1,'sum_fixedCost'} + genSum{1,'sum_variableCost'} - avgLMPtoLoad * genSum{1,'sum_annualGen'};
            envirSurp = -sum(genSum{1, {'sum_damCO2', 'sum_damNOx', 'sum_damSO2'}});             
            merchanSurplus = genSum{1,'sum_annualGen'} * (avgLMPtoLoad - avgLMPtoGen);
            totalSurplus = sum([consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus]);

            % Output the summary table
            outputs = table(consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus, totalSurplus);
        end

        %% sumSurplusCpp: sum surplus for cpp  case
        function [outputs] = sumSurplusCpp(genSum, genTable, busTable, caseInfo, result, yearInfo, iYear, tableGroup)
            % Calculate surplus
            % original obj value is hourly base

            % Map the result table to variable
            avgLMPtoLoad = genSum{1, 'avgLMPtoLoad'};
            avgLMPtoGen = genSum{1, 'avgLMPtoGen'};

            % Get the fuel results
            [fuelResults] = E4STResult.sumByfuel(genTable, busTable);

            % Get new wind and new solar
            idxNewWind = strcmp(result.fuel, 'wind') & result.newgen == 1;
            idxNewSolar = strcmp(result.fuel, 'solar') & result.newgen == 1;
            % Get the RE cost to build
            newWindCapCost = result.fixedcost(idxNewWind) - 5.542237443; % subtract cost to keep
            newSolarCapCost = result.fixedcost(idxNewSolar) - 3.845890411; % subtract cost to keep
            newWindCap = result.reserve.qty.Rp_pos(idxNewWind);
            newSolarCap = result.reserve.qty.Rp_pos(idxNewSolar);

            % Calculate CO2 payments based on CO2 emission price or cost of cap-trade program
            pvSubsidyFac = yearInfo.scalingSubsidyPV(iYear);
            windSubsidyFac = yearInfo.scalingSubsidyWind(iYear);
            reSubsidy = ((newSolarCapCost * pvSubsidyFac / (1 - pvSubsidyFac))' * newSolarCap + ...
                        (newWindCapCost * windSubsidyFac / (1 - windSubsidyFac))' * newWindCap) * 8760;

            multipliers = E4STResult.getMultipliersCPP(result);
            co2PayToGov = tableGroup{:, 'multiplier'}' * tableGroup{:, 'totalCO2'};
            co2PayFromProducer = co2PayToGov;
            % co2PayToGov = yearInfo.emissionPrice(iYear, 1) * genSum{:, 'sum_CO2'};
            % co2PayFromProducer = 0; % already included in the gencost 

            govRevenue = genSum{1,'sum_tax'} + co2PayToGov - reSubsidy;
            welProProfit = avgLMPtoGen * genSum{1,'sum_annualGen'} - genSum{1,'sum_fixedCost'} - genSum{1,'sum_variableCost'}...
                            - co2PayFromProducer;

            consumerSurp = -result.f * sum(caseInfo.hours) + genSum{1,'sum_fixedCost'} + genSum{1,'sum_variableCost'} - avgLMPtoLoad * genSum{1,'sum_annualGen'};
            envirSurp = -sum(genSum{1, {'sum_damCO2', 'sum_damNOx', 'sum_damSO2'}});             
            merchanSurplus = genSum{1,'sum_annualGen'} * (avgLMPtoLoad - avgLMPtoGen);
            totalSurplus = sum([consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus]);

            % Output the summary table
            outputs = table(consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus, totalSurplus);
        end      

        %% sumSurplusFastsim: sum surplus for fastsim  case
        function [outputs] = sumSurplusFastsim(genSum, genTable, busTable, caseInfo, result, yearInfo, iYear)
            % Calculate surplus
            % original obj value is hourly base

            % Map the result table to variable
            avgLMPtoLoad = genSum{1, 'avgLMPtoLoad'};
            avgLMPtoGen = genSum{1, 'avgLMPtoGen'};

            % Get the fuel results
            [fuelResults] = E4STResult.sumByfuel(genTable, busTable);

            % Get new wind and new solar
            idxNewWind = strcmp(result.fuel, 'wind') & result.newgen == 1;
            idxNewSolar = strcmp(result.fuel, 'solar') & result.newgen == 1;
            newWindCapCost = result.fixedcost(idxNewWind) - 5.542237443; % subtract cost to keep
            newSolarCapCost = result.fixedcost(idxNewSolar) - 3.845890411; % subtract cost to keep
            newWindCap = result.reserve.qty.Rp_pos(idxNewWind);
            newSolarCap = result.reserve.qty.Rp_pos(idxNewSolar);

            % Calculate CO2 payments based on CO2 emission price or cost of cap-trade program
            pvSubsidyFac = yearInfo.scalingSubsidyPV(iYear);
            windSubsidyFac = yearInfo.scalingSubsidyWind(iYear);
            reSubsidy = ((newSolarCapCost * pvSubsidyFac / (1 - pvSubsidyFac))' * newSolarCap + ...
                        (newWindCapCost * windSubsidyFac / (1 - windSubsidyFac))' * newWindCap) * 8760;
            co2PayToGov = yearInfo.emissionPrice(iYear, 1) * genSum{:, 'sum_CO2'};
            co2PayFromProducer = 0; % already included in the gencost 

            govRevenue = genSum{1,'sum_tax'} + co2PayToGov - reSubsidy;
            welProProfit = avgLMPtoGen * genSum{1,'sum_annualGen'} - genSum{1,'sum_fixedCost'} - genSum{1,'sum_variableCost'}...
                            - co2PayFromProducer;

            consumerSurp = -result.f * sum(caseInfo.hours) + genSum{1,'sum_fixedCost'} + genSum{1,'sum_variableCost'} - avgLMPtoLoad * genSum{1,'sum_annualGen'};
            envirSurp = -sum(genSum{1, {'sum_damCO2', 'sum_damNOx', 'sum_damSO2'}});             
            merchanSurplus = genSum{1,'sum_annualGen'} * (avgLMPtoLoad - avgLMPtoGen);
            totalSurplus = sum([consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus]);

            % Output the summary table
            outputs = table(consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus, totalSurplus);
        end  

         %% sumSurplusRtp: sum surplus for rtp  case
        function [outputs] = sumSurplusRtp(genSum, genTable, busTable, caseInfo, result, yearInfo, iYear)
            % Calculate surplus
            % original obj value is hourly base

            % Map the result table to variable
            avgLMPtoLoad = genSum{1, 'avgLMPtoLoad'};
            avgLMPtoGen = genSum{1, 'avgLMPtoGen'};

            % Get the fuel results
            [fuelResults] = E4STResult.sumByfuel(genTable, busTable);

            % Get new wind and new solar
            idxNewWind = strcmp(result.fuel, 'wind') & result.newgen == 1;
            idxNewSolar = strcmp(result.fuel, 'solar') & result.newgen == 1;
            newWindCapCost = result.fixedcost(idxNewWind) - 5.542237443; % subtract cost to keep
            newSolarCapCost = result.fixedcost(idxNewSolar) - 3.845890411; % subtract cost to keep
            newWindCap = result.reserve.qty.Rp_pos(idxNewWind);
            newSolarCap = result.reserve.qty.Rp_pos(idxNewSolar);

            % Calculate CO2 payments based on CO2 emission price or cost of cap-trade program
            pvSubsidyFac = yearInfo.emissionPrice(iYear, 1) * genSum{:, 'sum_CO2'};
            co2PayFromProducer = 0; % already included in the gencost 

            govRevenue = genSum{1,'sum_tax'} + co2PayToGov - reSubsidy;
            welProProfit = avgLMPtoGen * genSum{1,'sum_annualGen'} - genSum{1,'sum_fixedCost'} - genSum{1,'sum_variableCost'}...
                            - co2PayFromProducer;

            consumerSurp = -result.f * sum(caseInfo.hours) + genSum{1,'sum_fixedCost'} + genSum{1,'sum_variableCost'} - avgLMPtoLoad * genSum{1,'sum_annualGen'};
            envirSurp = -sum(genSum{1, {'sum_damCO2', 'sum_damNOx', 'sum_damSO2'}});             
            merchanSurplus = genSum{1,'sum_annualGen'} * (avgLMPtoLoad - avgLMPtoGen);
            totalSurplus = sum([consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus]);

            % Output the summary table
            outputs = table(consumerSurp, envirSurp, welProProfit, govRevenue, merchanSurplus, totalSurplus);
        end  

        %% getAllTable: get both genTable and busTable
        function [genRes, busRes] = getAllTable(mpc, offer, result, caseInfo, yearInfo, year)
            % Get the basic data table
            genRes = getGenRes(mpc, offer, result, caseInfo, yearInfo, year);
            busRes = getBusRes(result, caseInfo);        
        end

        % weightPrice: calculate load or generator weighted LMP; weighted by first column
        % first column: region; second column: load or generation; third column: price
        function [outputs] = weightPrice(arg)
            arg = [arg table(arg{:,2} .* arg{:,3})];
            funcall = @sum;
            sumArg = varfun(funcall, arg, 'GroupingVariables', 1);
            outputs = sumArg{:, 4} ./ sumArg{:, 2};
        end

        % calcLineprofit: calculate the profit of the line, the weighted
        % flow and positive flow (Quebec to NYC)
        function [profit, weightedflow, posFlow] = calcLineprofit(result)
            chline = E4STResult.getCHline(result);
            if length(chline) > 1
                profit = abs((chline(:,1) .* E4STResult.caseInfo.hours)') * abs(chline(:,2)); 
                weightedflow = abs(chline(:,1))' * E4STResult.caseInfo.afProb;
                pos_idx = chline(:,1) > 0;
                posFlow = chline(pos_idx, 1)' * E4STResult.caseInfo.afProb(pos_idx);
            else
                profit = 0;
                weightedflow = 0;
                posFlow = 0;
            end
        end     

        % getCHline: get results of the CH line
        function [outputs] = getCHline(result)
            outputs = zeros(E4STResult.caseInfo.nHours, 2);
            idx1 = find(result.base.bus(:,1) == 147829); % for NYC
            idx2 = find(result.base.bus(:,1) == 180708); % for quebec
            for i = 1: E4STResult.caseInfo.nHours
                if i == 1
                    dcline = result.base.dcline;
                    bus = result.base.bus;
                else
                    dcline = result.cont(i - 1).dcline;
                    bus = result.cont(i - 1).bus;
                end 
                % Check if ch line exists: 1-power flow; 2-lmp difference
                if dcline(end, 1) == 180708 && dcline(end, 2) == 147829
                    outputs(i, 1) = dcline(end, 4);
                    outputs(i, 2) = (bus(idx1, 14) - bus(idx2, 14)) / E4STResult.caseInfo.afProb(i);
                else
                    outputs = 0;
                    break;
                end
            end
        end

    end   
end

