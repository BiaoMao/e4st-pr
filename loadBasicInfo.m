function [basicInfo] = loadBasicInfo(basicInfoFile, verbose)
%% loadBasicInfo: load the basicinfo data to a struct data
%
%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

	% Set default argin
	if nargin < 2
		verbose = 1; % show a little debug information
	end

    infoTable = readtable(basicInfoFile);
    % Load each column into a data and remove NaN element
    for i = 1 : size(infoTable, 2)
        basicInfo.(infoTable.Properties.VariableNames{i}) = infoTable{:, i};
        idx = isnan(basicInfo.(infoTable.Properties.VariableNames{i}));
        basicInfo.(infoTable.Properties.VariableNames{i}) = ...
            basicInfo.(infoTable.Properties.VariableNames{i})(~idx, 1);
    end

    % Debug information
    if verbose == 1
        fprintf('%d columns of basic input data are loaded\n', size(infoTable, 2));
    end