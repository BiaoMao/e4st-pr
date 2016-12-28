function gencost = makeDrStep(mpc, verbose)
%% makeDrStep: make the step gencost for the demand response
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

	% Get the parameters
    elsty=caseInfo.drElasticity(iYear);
    disCost=caseInfo.drDiscost;
    priceVec=caseInfo.drPricevector;

    % Calculate the pivot point first
    pivotPoint = [loadInfo(:, 1) pivotPrice];
    [pLoads, pPrice] = DRBuilder.getDlPrices(mpc, pivotPoint);
    [defaultLoads, defaultPrices]=DRBuilder.getDlPrices(mpc,loadInfo);    

    % Calculate the starting point: [loads, basePrices] 
    loads = pLoads .* ((defaultPrices + disCost) ./ (pPrice + disCost)).^elsty;
            
    loads=-loads; % make negative for dl
    iDl=find(strcmp(mpc.genfuel(:,1),'dl')); % Find loads needed to make DR
    nDl=length(iDl); 


    %A=loads./(basePrices+disCost).^elsty;

    % Now we want to transform our input vector of price percentages into actual prices
    % for the last two steps, we want to set very high prices, 
    % unless the percentages are actually greater than these high prices ($5,000 and $10,000)
    stepPrices=zeros(nDl,10);
    % for i=1:10
    %     stepPrices(:,i)=(priceVec(i).*(defaultPrices+disCost))-disCost; 
    % end
    % stepPrices(:,9) = max(5000, stepPrices(:,9));
    % stepPrices(:,10) = max(15000, stepPrices(:,10)); 
    % stepPrices(:, 1) = min(pPrice, stepPrices(:,1)); % prevent two high starting price
    for i=1:10
        stepPrices(:,i) = priceVec(i).*defaultPrices; 
    end
    stepPrices(:,9) = max(5000, stepPrices(:,9));
    stepPrices(:,10) = max(15000, stepPrices(:,10)); 
    stepPrices(:, 1) = min(pPrice, stepPrices(:,1)); % prevent too high starting price


    
    %stepPrices(:,8) = min(stepPrices(:,9), stepPrices(:,8)); % Make it in decreasing order

    % For the vertexPower, we calculate the amount of power demanded 
    % at each price in the stepPrices
    % We set the power demanded at the top price level equal to 0
    % We need to calculate the actual vertex price and power
    vertexPower=zeros(nDl,10);
    vertexPrice=stepPrices;
    for i=1:9
        vertexPrice(:,i)=mean([stepPrices(:,i),stepPrices(:,i+1)],2);
        vertexPower(:,i) = loads.*((stepPrices(:,i)+disCost) ./ (defaultPrices+disCost)).^elsty;
        
        if ~isreal(vertexPower(:,i))
            %display(strcat('Huge Negative price'))
        end
    end
    vertexPower(:,10)=0;

    % Now we need to calculate the total costs, or f(p) in the gencost notation
    % We need to do this iteratively, as each (Power*Price) block must be added 
    % to the result of the previous TotalCostVector block
    totalCostVec=zeros(nDl,10);
    for i=9:-1:1
        totalCostVec(:,i)=totalCostVec(:,i+1)+(vertexPower(:,i)-vertexPower(:,i+1)).*...    
        vertexPrice(:,i);
%                 totalCostVec(:,i) = totalCostVec(:,i+1) + (vertexPower(:,i)-vertexPower(:,i+1)) .*...
%                  vertexPrice(:,i + 1);
        %(mean([stepPrices(:,i),stepPrices(:,i+1)],2));                                            
    end

    % Finally, we export the resultant gencost
    % this is p(MW),f(p) ($/hr)
    gencost=zeros(nDl,24);
    head=[1 0 0 10];
    gencost(:,1:4)=repmat(head,nDl,1);
    for i=1:10
        gencost(:,2*i+3)=vertexPower(:,i);
        gencost(:,2*i+4)=totalCostVec(:,i);
    end

    % Filter out gencost with negative LMP and zero loads
    % Get index of negative and zero prices and zero loads
    idxNegLoad = defaultPrices <= 0 | loads == 0 | pPrice<=0;
    idxNegGencost = iDl(idxNegLoad); 
    if any(idxNegLoad)
      gencost(idxNegLoad, :) = mpc.gencost(idxNegGencost, :); % We need to use different index
    end
    % Test if the gencost is good
    if ~isreal(gencost)
        display('Bad gencost')
        pause();
    end
end