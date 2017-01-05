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
        function [mpc] = setHydroCap0(mpc, caseInfo)
            define_constants;
            hydroMap = caseInfo.hydroMap{:, :};
            groups = size(hydroMap, 2);            
            n_hydro = length(find(hydroMap));            
            map = zeros(n_hydro, mpc.ng);
            cap = zeros(n_hydro, 1);
            coeff = ones(mpc.ng, 1);
            coeff_type = ones(n_hydro, 1);
            idx = 1;
            % Scale for each constraint groups
            for i = 1:groups
                idx_members = find(hydroMap(:,i));
                n_members = length(idx_members);
                
                % Scale the hydro capacity to real data
                mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * caseInfo.hydroScaling(i);

                % Set the low hydro for QC&NL
                if i == 4
                    curCap = sum(mpc.gen(idx_members, PMAX));
                    scaler = (curCap + caseInfo.deltaQCNL) / curCap;
                    mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                end

                % Set the low hydro for Ontario
                if i == 2
                    curCap = sum(mpc.gen(idx_members, PMAX));
                    scaler = (curCap + caseInfo.deltaON) / curCap;
                    mpc.gen(idx_members, PMAX) = mpc.gen(idx_members, PMAX) * scaler;
                end
                
                % Get the hydro CF
                cap(idx: idx+n_members-1) = mpc.gen(idx_members, PMAX) * caseInfo.hydroCf(i);
                % Re-build map for each hydro
                for j = 1:n_members
                    map(idx, idx_members(j)) = 1;
                    idx = idx + 1;
                end
            end
            mpc.total_output.map = map;
            mpc.total_output.cap = cap;
            mpc.total_output.coeff = coeff;
            mpc.total_output.type = coeff_type;
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