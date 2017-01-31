classdef HydroEI
%   Apply changes to Hydro in EI

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
        %% constrainHydro: set the hydro constraint cap
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

            % Extract the bus info
            define_constants;
            genBus = array2table(mpc.gen(:, 1), 'VariableNames', {'bus'});
            genBus = join(genBus, caseInfo.locationInfo);

            % For all US hydro
            hydroMap = strcmp(genBus{:, 'Nation'}, 'US') & strcmp(mpc.genfuel, 'hydro');

            % Ontario
            hydroMap = [hydroMap strcmp(genBus{:, 'State'}, 'ontario') & strcmp(mpc.genfuel, 'hydro')];

            % Manitoba
            hydroMap = [hydroMap strcmp(genBus{:, 'State'}, 'manitoba') & strcmp(mpc.genfuel, 'hydro')];
            
            % Quebec & NL
            hydroMap = [hydroMap (strcmp(genBus{:, 'State'}, 'quebec') | strcmp(genBus{:, 'State'}, 'newfoundland and labrador'))...
                     & strcmp(mpc.genfuel, 'hydro')];

            % New Brunswick & Nova Scotia
            hydroMap = [hydroMap (strcmp(genBus{:, 'State'}, 'new brunswick') | strcmp(genBus{:, 'State'}, 'nova scotia'))...
                     & strcmp(mpc.genfuel, 'hydro')];

            % SK
            hydroMap = [hydroMap (strcmp(genBus{:, 'State'}, 'saskatchewan'))...
                     & strcmp(mpc.genfuel, 'hydro')];

            groups = size(hydroMap, 2); % find all hydro groups 
            n_hydro = length(find(hydroMap)); % find all hydro in the map            
            map = zeros(n_hydro, size(mpc.gen, 1));
            cap = zeros(n_hydro, 1);
            coeff = ones(size(mpc.gen, 1), 1);
            coeff_type = ones(n_hydro, 1);
            idx = 1;

            % Scale for each constraint groups
            for i = 1:groups
                idx_members = find(hydroMap(:,i));
                n_members = length(idx_members);
                
                % % Scale the hydro capacity to real data
                % mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * caseInfo.hydroInfo{i, 'scaling'};

                % % Set the low hydro for QC&NL
                % if i == 4
                %     curCap = sum(mpc.gen(idx_members, PMAX));
                %     scaler = (curCap + caseInfo.deltaQCNL) / curCap;
                %     mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                % end

                % % Set the low hydro for Ontario
                % if i == 2
                %     curCap = sum(mpc.gen(idx_members, PMAX));
                %     scaler = (curCap + caseInfo.deltaON) / curCap;
                %     mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                % end
                
                % Get the hydro CF
                cap(idx: idx+n_members-1) = mpc.gen(idx_members, PMAX) * caseInfo.hydroInfo{i, 'Cf'};
                % Re-build map for each hydro
                for j = 1:n_members
                    map(idx, idx_members(j)) = 1;
                    idx = idx + 1;
                end

                % Set hydro AFs
                mpc.availability_factor(idx_members, :) = caseInfo.hydroInfo{i, 'Af'};

                % Set Pmin of Hydro
                mpc.gen(idx_members, PMIN) = mpc.gen(idx_members, PMAX) * caseInfo.hydroInfo{i, 'Pmin'};
           
                % Debug information
                if verbose == 1
                    fprintf('%d hydro constraints in Gourp %d are set in EI\n', n_members, i);
                end
            end
            mpc.total_output.map = [mpc.total_output.map; map];
            mpc.total_output.cap = [mpc.total_output.cap; cap];
            mpc.total_output.coeff(:, 1) = coeff; % first column is all ones
            mpc.total_output.type = [mpc.total_output.type; coeff_type];  
        end

        %% setHydroCap10: definitely correct hydro cap in CA in 2025
        function [mpc, offer] = setHydroCap10(mpc, offer, caseInfo)
            define_constants;
            groups = size(hydroMap, 2);
            deltaQCNL = 1998;
            deltaMan = 695;
            deltaSK = 50;
 
            % Scale for each constraint groups
            for i = 1:groups
                idx_members = find(hydroMap(:,i));

                % Set the low hydro for QC&NL
                if i == 4
                    curCap = sum(mpc.gen(idx_members, PMAX));
                    scaler = (curCap + deltaQCNL) / curCap;
                    mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                    offer(idx_members, 2) = offer(idx_members, 2) * scaler;
                end

                % Set the low hydro for Manitoba
                if i == 3
                    curCap = sum(mpc.gen(idx_members, PMAX));
                    scaler = (curCap + deltaMan) / curCap;
                    mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                    offer(idx_members, 2) = offer(idx_members, 2) * scaler;
                end                

                % Set the low hydro for Sk
                if i == 6
                    curCap = sum(mpc.gen(idx_members, PMAX));
                    scaler = (curCap + deltaSK) / curCap;
                    mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                    offer(idx_members, 2) = offer(idx_members, 2) * scaler;
                end 
            end
        end

        %% setHydroCap10High: correct hydro cap in CA in 2025 for high hydro case
        function [mpc, offer] = setHydroCap10High(mpc, offer, caseInfo)
            define_constants;
            groups = size(hydroMap, 2);
            deltaQCNL = 2264;
            deltaOntario = 25;
 
            % Scale for each constraint groups
            for i = 1:groups
                idx_members = find(hydroMap(:,i));

                % Set the low hydro for QC&NL
                if i == 4
                    curCap = sum(mpc.gen(idx_members, PMAX));
                    scaler = (curCap + deltaQCNL) / curCap;
                    mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                    offer(idx_members, 2) = offer(idx_members, 2) * scaler;
                end

                % Set the low hydro for Ontario
                if i == 2
                    curCap = sum(mpc.gen(idx_members, PMAX));
                    scaler = (curCap + deltaOntario) / curCap;
                    mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                    offer(idx_members, 2) = offer(idx_members, 2) * scaler;
                end           
            
            end
        end

        %% setHydroCap20High: correct hydro cap in CA in 2035 for high hydro case
        function [mpc, offer] = setHydroCap20High(mpc, offer, caseInfo)
            define_constants;
            groups = size(hydroMap, 2);
            deltaQCNL = 1200;
            deltaMan = 1485;
 
            % Scale for each constraint groups
            for i = 1:groups
                idx_members = find(hydroMap(:,i));

                % Set the low hydro for QC&NL
                if i == 4
                    curCap = sum(mpc.gen(idx_members, PMAX));
                    scaler = (curCap + deltaQCNL) / curCap;
                    mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                    offer(idx_members, 2) = offer(idx_members, 2) * scaler;
                end

                % Set the low hydro for Manitoba
                if i == 3
                    curCap = sum(mpc.gen(idx_members, PMAX));
                    scaler = (curCap + deltaMan) / curCap;
                    mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                    offer(idx_members, 2) = offer(idx_members, 2) * scaler;
                end           
            
            end
        end     

        %% configHydro: re-write the availability factors and minimum dispatch
        function [mpc, offer] = setHydroAF(mpc, offer, caseInfo)
            groups = size(caseInfo.hydroAf, 2);
            for i = 1:groups
                idx_members = find(hydroMap(:,i));
                mpc.availability_factor(idx_members, :) = caseInfo.hydroAf(i);
                offer(idx_members, 13) = caseInfo.hydroPmin(i) / caseInfo.hydroAf(i);
            end
        end

       
    end    
end