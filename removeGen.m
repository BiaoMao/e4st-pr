function [mpc, offer] = removeGen(mpc, offer, idx)
%% removeGen: remove a set of gen from mpc and offer
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    mpc.gen(idx, :) = [];
    mpc.gencost(idx, :) = [];
    mpc.genfuel(idx, :) = [];
    mpc.gen_aux_data(idx, :) = [];
    mpc.newgen(idx, :) = [];
    mpc.availability_factor(idx, :) = [];
    offer(idx, :) = [];
    mpc.ng = size(mpc.gen, 1); 

    % Update total_output constraints
    if isfield(mpc, 'total_output') && ~isempty(mpc.total_output)
        mpc.total_output.map(:, idx) = [];
        mpc.total_output.coeff(idx, :) = [];
    end

    % Update capacity constraints
    if isfield(mpc, 'caplim') && ~isempty(mpc.caplim)
        mpc.caplim.map(:, idx) = [];
    end
end