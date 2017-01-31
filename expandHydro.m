function mpc = expandHydro(mpc, caseInfo, verbose)
%Expand hydro capacity based on States

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
    genBus = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
    genBus = join(genBus, caseInfo.locationInfo);
    % Update hydro cap by State
    nStates = size(caseInfo.hydroAdd, 1);
    for i = 1 : nStates        
        idx = strcmp(mpc.genfuel, 'hydro') & strcmp(genBus{:, 'State'}, caseInfo.hydroAdd{i, 'State'});
        curCap = sum(mpc.gen(idx, PMAX));
        scaler = (curCap + caseInfo.hydroAdd{i, 'capAdd'}) / curCap;
        mpc.gen(idx, PMAX) = mpc.gen(idx, PMAX) * scaler;
        % Debug information
        if verbose == 1
            fprintf('Hydro in %s is expanded by %f\n', caseInfo.hydroAdd{i, 'State'}{:}, scaler);
        end
    end
end





 