function [caseInfo] = scaleCost2Build(mpc, caseInfo, yearInfo, year, verbose)
%  Apply scalings to cost to build

%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

    % Set default argin
    if nargin < 5
        verbose = 1; % show a little debug information
    end

    % Input data
    idxYear = find(strcmp(yearInfo.Properties.VariableNames, ['Y', num2str(year)]));

    % Scale cost to build for solar
    if any(strcmp(yearInfo.Properties.RowNames, 'c2bSolar'))
        scaling = yearInfo{'c2bSolar', idxYear} / yearInfo{'c2bSolar', idxYear - 1};
        caseInfo.genInfo{'solar', 'Cost2Build'} = caseInfo.genInfo{'solar', 'Cost2Build'} * scaling;
        if verbose == 1
            fprintf('Scale cost to build for solar by %f\n', yearInfo{'c2bSolar', idxYear});
        end
    end

    % Scale cost to build for wind
    if any(strcmp(yearInfo.Properties.RowNames, 'c2bWind'))
        scaling = yearInfo{'c2bWind', idxYear} / yearInfo{'c2bWind', idxYear - 1};
        caseInfo.genInfo{'wind', 'Cost2Build'} = caseInfo.genInfo{'wind', 'Cost2Build'} * scaling;
        if verbose == 1
            fprintf('Scale cost to build for wind by %f\n', yearInfo{'c2bWind', idxYear});
        end
    end

    % Scalie cost to build for other types
    if any(strcmp(yearInfo.Properties.RowNames, 'c2bOther'))
        scaling = yearInfo{'c2bOther', idxYear} / yearInfo{'c2bOther', idxYear - 1};
        idxOther = ~(strcmp(caseInfo.genInfo.Properties.RowNames, 'solar') | ...
                    strcmp(caseInfo.genInfo.Properties.RowNames, 'wind'));
        caseInfo.genInfo{idxOther, 'Cost2Build'} = caseInfo.genInfo{idxOther, 'Cost2Build'} * scaling;
        if verbose == 1
            fprintf('Scale cost to build for others by %f\n', yearInfo{'c2bOther', idxYear});
        end
    end
