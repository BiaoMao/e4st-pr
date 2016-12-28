function [offer] = makeGenOffer(mpc, genInfo, verbose)
%% makeGenOffer: set offer table for all existing generators
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

    % Set default argin
    if nargin < 3
    	verbose = 1; % show a little debug information
    end

    define_constants;
    % Intialize the offer table
    offer = zeros(mpc.ng, 13);
    offer(:, 3) = 0;
    offer(:, 4) = Inf;

    % Offer table coefficient/scalling
    coff = 1;

    % Get all the fuelltypes
    fuels = unique(mpc.genfuel);

    for fuel = fuels'
        iGen = strcmp(mpc.genfuel, fuel);
        % Set the PositiveActiveReservePrice and PositiveActiveReserveQuantity
        % Add cost to keep, tax and insurance to fixed cost
        posPrice = sum(genInfo{fuel, {'Cost2Keep', 'Tax', 'Insurance'}}, 2);
        posQty = max(mpc.gen(iGen, PMAX).*coff, -mpc.gen(iGen, PMIN).*coff);
        offer(iGen, 1) = posPrice;
        offer(iGen, 2) = posQty;
    end

    % Debug information
    if verbose == 1
        fprintf('Set offer for all existing generators \n');
    end

 