function [mpc] = updateConsDim(mpc)
%% updateConsDim: update constraints dims after adding generators
%   Dims include total_output and caplim
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

    % Update dims of total_output
    if isfield(mpc, 'total_output')
        [n_contraints, old_ng] = size(mpc.total_output.map);
        delta_ng = mpc.ng - old_ng;
        % Add all zeros to map and coeff
        % Will rebuild total_output every time
        % Not necessary change here
        if delta_ng > 0
            mpc.total_output.map = [mpc.total_output.map zeros(n_contraints, delta_ng)];
            mpc.total_output.coeff = [mpc.total_output.coeff;...
                    zeros(delta_ng, size(mpc.total_output.coeff, 2))];
        end
    end

    % Update dims of caplim
    if isfield(mpc, 'caplim')
        [n_contraints, old_ng] = size(mpc.caplim.map);
        delta_ng = mpc.ng - old_ng;
        % Add all zeros to map
        if delta_ng > 0
            mpc.caplim.map = [mpc.caplim.map zeros(n_contraints, delta_ng)];
        end
    end
 