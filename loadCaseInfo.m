function [caseInfo] = loadCaseInfo(caseInfoFile, verbose)
%% loadBasicInfo: load the caseInfo data to a struct data
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

    infoTable = readtable(caseInfoFile);
    % Load each column into a data and remove NaN/empty element
    for i = 1 : size(infoTable, 2)
        structName = infoTable.Properties.VariableNames{i};
        vector = infoTable{:, i};
        if iscell(vector)
            idx = isempty(vector);   
            caseInfo.(structName) = vector(~idx);         
        else 
            idx = isnan(vector);
            caseInfo.(structName) = vector(~idx);
        end      
    end

    % Debug information
    if verbose == 1
        fprintf('%d columns of caseInfo data are loaded\n', size(infoTable, 2));
    end