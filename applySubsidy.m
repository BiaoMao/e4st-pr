function [mpc] = applySubsidy(mpc, yearInfo, year, verbose)
%% applySubsidy: Apply subsidies for certain types of generators
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

    % Set default argin
    if nargin < 4
    	verbose = 1; % show a little debug information
    end

    % Format the year as a string
	year = ['Y', num2str(year)];
	idxYear = find(strcmp(yearInfo.Properties.VariableNames, year));  

    % Apply wind subsidy
    fuelType = 'wind';
    if any(strcmp(yearInfo.Properties.RowNames, 'windSub'))
        % For non-new gen
        idx = strcmp(mpc.genfuel, fuelType) & mpc.newgen ~= 1;      
        % The subsidy change from last year
        delta = yearInfo{'windSub', idxYear} - yearInfo{'windSub', idxYear - 1};
        mpc.gencost(idx, 5) = mpc.gencost(idx, 5) - delta;  

        % For new gen
        idx = strcmp(mpc.genfuel, fuelType) & mpc.newgen == 1; 
        mpc.gencost(idx, 5) = mpc.gencost(idx, 5) - yearInfo{'windSub', idxYear}; 

        mpc.wind_subsidy = yearInfo{'windSub', idxYear};

        if verbose == 1
            fprintf('Apply $%f subsidy to gencost for %s\n', ...
                yearInfo{'windSub', idxYear}, fuelType);
        end 
    end