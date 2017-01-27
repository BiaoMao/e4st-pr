classdef HydroWECC
%   Apply changes to Hydro in WECC

%   E4ST
%   Copyright (c) 2012-2016 by Power System Engineering Research Center (PSERC)
%   by Biao Mao, Rensselaer Polytechnic Institute
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.
    
    properties  (Constant)

    end

    methods (Static)
        function [mpc] = setHydroAfCf(mpc, caseInfo, verbose)
            % Set default argin
            if nargin < 3
                verbose = 1; % show a little debug information
            end

            % Check if total_output filed is initialized
            if ~isfield(mpc, 'total_output') || isempty(mpc.total_output)
                mpc.total_output.map = [];
                mpc.total_output.cap = [];
                mpc.total_output.coeff = [];
                mpc.total_output.type = [];
            end

            define_constants;
            
            % Extract the bus info
            genBus = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
            genBus = join(genBus, caseInfo.locationInfo);

            % For all US hydro
            hydroMap = strcmp(genBus{:, 'Nation'}, 'US') & strcmp(mpc.genfuel, 'hydro');

            % British Columbia (BC)
            hydroMap = [hydroMap strcmp(genBus{:, 'State'}, 'british columbia') & strcmp(mpc.genfuel, 'hydro')];

            % Alberta (AB)
            hydroMap = [hydroMap strcmp(genBus{:, 'State'}, 'alberta') & strcmp(mpc.genfuel, 'hydro')];
            
            groups = size(hydroMap, 2); % find all hydro groups 
            n_hydro = length(find(hydroMap)); % find all hydro in the map            
            map = zeros(n_hydro, size(mpc.gen, 1));
            cap = zeros(n_hydro, 1);
            coeff = ones(size(mpc.gen, 1), 1);
            coeff_type = ones(n_hydro, 1);
            idx = 1;
            for i = 1 : groups
                idx_members = find(hydroMap(:,i)); % use the actual index instead of binary 
                n_members = length(idx_members);
                
                % Set hydro CFs
                cap(idx: idx+n_members-1) = mpc.gen(idx_members, PMAX) * caseInfo.hydroInfo{i, 'Cf'};
                % Re-build map for each hyro
                for j = 1:n_members
                    map(idx, idx_members(j)) = 1;
                    idx = idx + 1;
                end

                % Set hydro AFs
                mpc.availability_factor(idx_members, :) = caseInfo.hydroInfo{i, 'Af'};

                % Set Pmin of Hydro
                mpc.gen(idx_members, PMIN) = mpc.gen(idx_members, PMAX) * caseInfo.hydroInfo{i, 'Pmin'};
            end
            mpc.total_output.map = [mpc.total_output.map; map];
            mpc.total_output.cap = [mpc.total_output.cap; cap];
            mpc.total_output.coeff(:, 1) = coeff; % first column is all ones
            mpc.total_output.type = [mpc.total_output.type; coeff_type];            

            % Debug information
            if verbose == 1
                fprintf('Hydro AFs, Pmins and CFs are set\n');
            end
        end

        function mpc = expandHydro(mpc, caseInfo, verbose)
            % Set default argin
            if nargin < 3
                verbose = 1; % show a little debug information
            end

            % Update hydro cap in Washinton state
            genBus = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
            genBus = join(genBus, caseInfo.locationInfo);
            idxWA = strcmp(mpc.genfuel, 'hydro') & strcmp(genBus{:, 'State'}, 'washinton');
            curCap = sum(mpc.gen(idxWA, PMAX));
            scaler = (curCap + caseInfo.deltaHydroWA) / curCap;
            mpc.gen(idxWA, PMAX) = mpc.gen(idxWA, PMAX) * scaler;
             % Debug information
            if verbose == 1
                fprintf('Hydro in Washinton is expanded by %f MW\n', scaler);
            end

            % Update hydro cap in BC   
            idxWA = strcmp(mpc.genfuel, 'hydro') & strcmp(genBus{:, 'State'}, 'british columbia');
            curCap = sum(mpc.gen(idxWA, PMAX));
            scaler = (caseInfo.hydroBC) / curCap;
            mpc.gen(idxWA, PMAX) = mpc.gen(idxWA, PMAX) * scaler;
             % Debug information
            if verbose == 1
                fprintf('Hydro in british columbia is expanded by %f MW\n', scaler);
            end
        end

    end
end



 