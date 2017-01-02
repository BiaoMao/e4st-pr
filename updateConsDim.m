function [mpc] = updateConsDim(mpc)
%% updateConsDim: updated constraints dimentions after adding generators
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

    [n_contraints, old_ng] = size(mpc.total_output.map);
    delta_ng = mpc.ng - old_ng;
    mpc.total_output.map = [mpc.total_output.map zeros(n_contraints, delta_ng)];        
    
    mpc.total_output.coeff = [mpc.total_output.coeff; ones(delta_ng, 1)];

 