function [caseInfo] = applyPVSubsidy(caseInfo, yearInfo, year, verbose)
%  Apply scalings to cost to build for solar subsidy

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
    if any(strcmp(yearInfo.Properties.RowNames, 'solarSub'))
        if idxYear == 1
            scaling = yearInfo{'solarSub', idxYear};
        else 
            scaling = yearInfo{'solarSub', idxYear} / yearInfo{'solarSub', idxYear - 1};
        end
        caseInfo.genInfo{'solar', 'Cost2Build'} = caseInfo.genInfo{'solar', 'Cost2Build'} * scaling;
        if verbose == 1
            fprintf('Scale cost to build for solar by %f for subsidy\n', yearInfo{'solarSub', idxYear});
        end
    end
