function saveRawRes(fileName, mpc, offer, result, caseInfo, yearInfo)
%% updateDamage: update locational damages for NOx and SO2
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

    curYear = ['Y' num2str(result.year)];
    fVal = result.(curYear).opf_results.f;
    result.(curYear).opf_results = [];
    result.(curYear).opf_results.f = fVal;
    save(fileName, 'mpc', 'offer', 'result', 'caseInfo', 'yearInfo', '-v7.3');
end