function contab = makeContDr(mpc, caseInfo, yearInfo, busRes, year, mode, verbose)
%% makeContDr: make the step gencost of the demand response for contingency hours
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    % Set default argin
    if nargin < 7
        verbose = 1; % show a little debug information
    end

	nHours = caseInfo.nHours; %total hours(including base hours)
    iDl = find(strcmp(mpc.genfuel(:, 1), 'dl'));
    nDl = length(iDl);

    % The rows generated by 'dl' is nDl * 18 * (nHours - 1); 
    % 10 steps of pairs; last two are zero
    nConts = nDl * 18 * (nHours - 1); % total rows in contab
    contab = zeros(nConts, 7);

    % Save the original mpc for scalling
    for i = 2 : nHours            
        gencost = makeDrStep(mpc, caseInfo, yearInfo, busRes, i, year, mode);
        iCont = i - 1;

        % Create contab for one hour
        col2Change = 1 : 22; % colunms need to be changed
        nCol = length(col2Change);% total columns
        nSigCont = nDl * nCol;% total rows of contabs for one hour
        gencostVec = reshape(gencost(:, col2Change),[],1); % value need to be changed
        headVec1 = repmat([iCont caseInfo.probability(i) 9],nSigCont,1);
        headVec2 = repmat(iDl,nCol,1);
        headVec3 = reshape(repmat(col2Change,nDl,1),[],1); % column id to change
        headVec4 = repmat([1],nSigCont,1);
        contab1 = [headVec1 headVec2 headVec3 headVec4 gencostVec];
        contab(1+(iCont-1)*nSigCont:iCont*nSigCont,:) = contab1;

        % Add load limits using the first step load
        % iUpdate = find(gencost(: ,5) < 0);
        % nUpdate = length(iUpdate);
        % headVec1 = repmat([iCont caseInfo.probability(i) 2],nUpdate,1);
        % headVec2 = iDl(iUpdate);
        % headVec3 = repmat([10 1], nUpdate, 1); % column id to change, 10: load limits, 1: new values
        % headVec4 = gencost(iUpdate, 5); % the firs step of loads

        % contab = [contab; [headVec1 headVec2 headVec3 headVec4]];
    end
end