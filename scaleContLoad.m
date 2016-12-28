function [contab] = scaleContLoad(mpc, basicInfo, scalings, verbose)
%% scaleContLoad: scale loads in a area with scallings in contingency hour
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

    % Initialize 
    define_constants;
    CT_LOAD = 4; % modify fixed and dispatchable loads
    CT_TABLE = 8; % area-wide load changes
    nHours = length(basicInfo.hours);
    area = unique(mpc.bus(:, BUS_AREA));
    conNum = length(area) * (nHours-1);
    contab = zeros(conNum, 7);
    i = 1;
    for iHour = 2 : nHours  
      for iArea = 1 : length(area)
          contab(i,:) = [iHour - 1 basicInfo.probability(iHour) CT_TABLE area(iArea) CT_LOAD  2 ...
           scalings{iArea, iHour}];
          i = i + 1;
      end
    end

    % Debug information
    if verbose == 1
        fprintf('Scale loads in %d area in %d contingency hour\n', length(area), nHours - 1);
    end

 