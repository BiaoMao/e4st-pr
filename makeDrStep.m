function gencost = makeDrStep(mpc, caseInfo, yearInfo, busRes, hour, year, mode, verbose)
%% makeDrStep: make the step gencost of the demand response
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    % Set default argin
    if nargin < 8
        verbose = 1; % show a little debug information
    end

	% Get the parameters
    curYear = ['Y' num2str(year)];
    elsty = yearInfo{'elasticity', curYear};
    disCost = caseInfo.discost;
    priceVec = caseInfo.priceVec;

    % Decide what price will be used as pivot price
    if strcmp(mode, 'annual')
        pPrice = busRes.sysLMP;
    end

    % Calculate the pivot points and defualt prices
    % Convert to gen bus index first
    busTable = array2table([mpc.bus(:, 1), busRes.loadHourly(:, hour), busRes.lmpHourly(:, hour)],...
        'VariableNames', {'bus', 'load', 'lmp'});
    genTable = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
    genTable = join(genTable, busTable);
    pLoads = genTable{:, 'load'};    
    defaultPrices = genTable{:, 'lmp'};
    
    % Select 'dl' only
    iDl = find(strcmp(mpc.genfuel(:, 1), 'dl')); % Find loads needed to add step gencost
    nDl = length(iDl); 
    pLoads = pLoads(iDl, 1);
    defaultPrices = defaultPrices(iDl, 1);

    % Calculate the starting point: [loads, basePrices] 
    loads = pLoads .* ((defaultPrices + disCost) ./ (pPrice + disCost)).^elsty;
            
    loads = -loads; % make negative for 'dl'

    % Now we want to transform our input vector of price percentages into actual prices
    % for the last two steps, we want to set very high prices, 
    % unless the percentages are actually greater than these high prices ($5,000 and $10,000)
    % use 10 steps of prices
    stepPrices = zeros(nDl, 10);
    for i = 1 : 10
        stepPrices(:,i) = priceVec(i) .* defaultPrices; 
    end
    stepPrices(:,9) = max(5000, stepPrices(:,9));
    stepPrices(:,10) = max(15000, stepPrices(:,10)); 
    stepPrices(:, 1) = min(pPrice, stepPrices(:,1)); % prevent too high starting price    

    % For the vertexPower, we calculate the amount of power demanded 
    % at each price in the stepPrices
    % We set the power demanded at the top price level equal to 0
    % We need to calculate the actual vertex price and power
    vertexPower = zeros(nDl,10);
    vertexPrice = stepPrices;
    for i = 1 : 9
        vertexPrice(:, i) = mean([stepPrices(:, i), stepPrices(:, i + 1)], 2);
        vertexPower(:, i) = loads .* ((stepPrices(:, i) + disCost) ./ (defaultPrices + disCost)) .^ elsty;        
        if ~isreal(vertexPower(:, i))
            fprintf('Huge negative price at hour %d\n', hour);
        end
    end
    vertexPower(:, 10) = 0;

    % Now we need to calculate the total costs, or f(p) in the gencost notation
    % We need to do this iteratively, as each (Power*Price) block must be added 
    % to the result of the previous TotalCostVector block
    totalCostVec = zeros(nDl, 10);
    for i = 9 : -1 : 1
        totalCostVec(:, i) = totalCostVec(:, i + 1) +...
            (vertexPower(:, i) - vertexPower(:, i + 1)) .* vertexPrice(:, i);                                          
    end

    % Finally, we export the resultant gencost
    % this is p(MW),f(p) ($/hr)
    gencost = zeros(nDl, 24);
    head=[1 0 0 10];
    gencost(:, 1 : 4) = repmat(head, nDl, 1);
    for i = 1 : 10
        gencost(:, 2 * i + 3) = vertexPower(:, i);
        gencost(:, 2 * i + 4) = totalCostVec(:, i);
    end

    % Filter out gencost with negative LMP and zero loads
    % Get index of negative and zero prices and zero loads
    idxNegLoad = defaultPrices <= 0 | loads == 0 | defaultPrices == 5000 | pLoads <= 0;
    idxNegGencost = iDl(idxNegLoad); 
    if any(idxNegLoad)
        % Recover to original cost
        gencost(idxNegLoad, 5) = mpc.gencost(idxNegGencost, 5); % We need to use different index
        gencost(idxNegLoad, 6 : end) = 0;
    end
    % Test if the gencost is good
    if ~isreal(gencost)
        fprintf('Bad gencost at hour %d\n', hour);
        pause();
    end
end