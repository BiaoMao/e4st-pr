function [result] = removeRes(result, idx)
%% removeRes: remove a set of gen from result
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    result.base.gen(idx, :) = [];
    nCont = size(result.cont, 2);
    for i = 1 : nCont
        result.cont(1, i).gen(idx, :) = [];
    end
    result.energy.sum_muPmax(idx, :) = [];
    result.energy.sum_muPmin(idx, :) = [];
    result.energy.Pc(idx, :) = [];
    result.energy.Gmax(idx, :) = [];
    result.energy.Gmin(idx, :) = [];
    result.energy.Qmax(idx, :) = [];
    result.energy.Qmin(idx, :) = [];
    result.energy.mu.alphaP(idx, :) = [];
    result.energy.mu.Pc(idx, :) = [];

    result.reserve.mu.Rp_pos(idx, :) = [];
    result.reserve.mu.Rp_neg(idx, :) = [];
    result.reserve.mu.Rpmax_pos(idx, :) = [];
    result.reserve.mu.Rpmax_neg(idx, :) = [];
    result.reserve.qty.Rp_pos(idx, :) = [];
    result.reserve.qty.Rp_neg(idx, :) = [];
    result.reserve.prc.Rp_pos(idx, :) = [];
    result.reserve.prc.Rp_neg(idx, :) = [];
end