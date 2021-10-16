classdef fluid
    %FLUID: Stores fluid properties on each cell.
    
    properties
        %% Properties at time level n+1
        vecP; % Oil pressure[psi] and water saturation [Po,Sw...]: 2 by 2*nCell
              % saves values for two time levels: n+1 at k iter., n+1
              % at k+1 iterations.
        Po; %Oil pressure[psi](= vecP(1:2:end)): 2 by nCell
        Sw; %Water saturation(= vecP(2:2:end)): 2 by nCell
        Bo; % Oil formation volume factor [dimensionless; reservor V / STD V]
        Bw; % Water formation volume factor [dimensionless; reservor V / STD V]
        bo; % 1/Oil formation volume factor [dimensionless; reservor V / STD V]
        bw; % 1/Water formation volume factor [dimensionless; reservor V / STD V]
        gamma_denO; % Oil density multiplied by g/gc [psi/ft]
        gamma_denW; % Water density multiplied by g/gc [psi/ft]
        visO; % Oil viscosity at reservoir condition [cp]
        visW; % Water viscosity at reservoir condition [cp]
        Kro % Relative permeability of oil
        Krw; % Relative permeability of water
        denO_std; % Oil density at standard condition [lbm/ft3]
        denW_std; % Water density at standard condition [lbm/ft3]
        visO_std; % Oil viscosity at standard condition [cp]
        visW_std; % Water viscosity at standard condition [cp]
        Pref = 14.7; % Reference pressure[psi]
        conv_beta = 1/144; %ft2/in2
        
        %% Properties at time level n
        vecPPrev; %Oil pressure[psi] and water saturation over time [Po,Sw...]
        PoPrev; %Previous oil pressure[psi]
        SwPrev; % Previous water saturation
        boPrev; % 1/Oil formation volume factor [dimensionless; reservor V / STD V]
        bwPrev; % 1/Water formation volume factor [dimensionless; reservor V / STD V]
        
        %% Initialization of properties at time n+1 at previouse iter.
        PoPrevIter;
        SwPrevIter;
        
        %% Derivatives of properties at time level n+1
        deriv_bo; % d(Bo)/d(Po) at time level n+1
        deriv_bw; % d(Bw)/d(Po) at time level n+1
        derivVisO; % d(visO)/d(Po) at time level n+1
        derivVisW; % d(visW)/d(Po) at time level n+1dx
        derivKro; % d(Kro)/d(Sw) at time level n+1
        derivKrw; % d(Krw)/d(Sw) at time level n+1
        compOil = 2e-5; % Oil compressibility [psi-1]
        compWat % Water compressibility [psi-1]        
        compVisOil; % Water compressibility [psi-1]
        compVisWat; % Water compressibility [psi-1]
    end
    
    methods
        %% Class initialization
        function obj = fluid(grid)
            nTot = grid.nCell;
            % Initialization of properties at time n+1
            obj.vecP = zeros(1, 2*nTot);
            obj.Po = zeros(1, nTot);
            obj.Sw = zeros(1, nTot);
            obj.bo = zeros(1, nTot);
            obj.bw = zeros(1, nTot);
            obj.visO = zeros(1, nTot);
            obj.visW = zeros(1, nTot);
            obj.Kro = zeros(1, nTot);
            obj.Krw = zeros(1, nTot);
            obj.gamma_denO = zeros(1, nTot);
            obj.gamma_denW = zeros(1, nTot);
            
            % Initialization of properties at time n
            obj.vecPPrev = zeros(1, 2*nTot);
            obj.PoPrev = zeros(1, nTot);
            obj.SwPrev = zeros(1, nTot);
            obj.boPrev = zeros(1, nTot);
            obj.bwPrev = zeros(1, nTot);
            
            % Initialization of properties at time n+1 at previouse iter.
            obj.PoPrevIter = zeros(1, nTot);
            obj.SwPrevIter = zeros(1, nTot);
            
            % Initialization of derivatives 
            obj.deriv_bo = zeros(1, nTot);
            obj.deriv_bw = zeros(1, nTot);
            obj.derivVisO = zeros(1, nTot);
            obj.derivVisW = zeros(1, nTot);
            obj.derivKro = zeros(1, nTot);
            obj.derivKrw = zeros(1, nTot);            
        end
        
        %% Initialize pressure using bisection method
        function [vecP, Po] = initializePressure(obj, grid)
            errLimit = 1e-16;
            Ptemp = obj.Po(1)*zeros(1,grid.nz); 
            for k = 1:grid.nz
                err = 1;
                if k == 1
                    Plow = obj.Po(1);
                else
                    Plow = Ptemp(k-1);
                end
                iCell = 1+(k-1)*grid.nx*grid.ny;
                Phigh = Plow + 3*obj.denO_std/144*grid.dz(iCell);
                while abs(err) > errLimit
                    Pestimated = (Phigh+Plow)/2;           
                    gamma_denO(k) = obj.calcOilDensity(Pestimated); 
                    if k == 1
                        Ptemp(k) = obj.Po(k) + gamma_denO(k)*grid.dz(iCell)/2;
                    else
                        Ptemp(k) = Ptemp(k-1)+(gamma_denO(k-1)+gamma_denO(k))/2*grid.dz(iCell);
                    end
                    err = Ptemp(k) - Pestimated; 
                    if err > 0
                        Plow = Pestimated;
                    else
                        Phigh = Pestimated;
                    end
                end                
                Po(iCell:(iCell+grid.nx*grid.ny-1)) = Ptemp(k);
            end
            vecP = obj.vecP;
            vecP(1:2:end) = Po;
        end
        
        %% Local function for initializing pressure.
        function gamma_denO = calcOilDensity(obj, p)
            bo = exp(obj.compOil*(p-obj.Pref));
            gamma_denO = obj.conv_beta*obj.denO_std.*bo;
        end
        
        %% Update [bw, bo, gamma_denO, gamma_denW]
        function [bw, bo, gamma_denW, gamma_denO, Bw, Bo] = updateFVF(obj)
            Bo = exp(-obj.compOil*(obj.Po-obj.Pref));
            Bw = exp(-obj.compWat*(obj.Po-obj.Pref));
            bw = 1./Bw;
            bo = 1./Bo;
            gamma_denW = obj.conv_beta*obj.denW_std.*bw;
            gamma_denO = obj.conv_beta*obj.denO_std.*bo; 
        end
        
        %% Update [bw, bo]
        function [bw, bo] = updatebwbo(obj, Po)
            Bo = exp(-obj.compOil*(Po-obj.Pref));
            Bw = exp(-obj.compWat*(Po-obj.Pref));
            bw = 1./Bw;
            bo = 1./Bo; 
        end
        
        %% Update viscosity[visW, visO]
        function [visW, visO] = updateVis(obj)
            visW = obj.visW_std*exp(obj.compVisWat*(obj.Po-obj.Pref));
            visO = obj.visO_std*exp(obj.compVisOil*(obj.Po-obj.Pref));
        end
        
        %% Update relative permeabilities[Krw, Kro]
        function [Krw, Kro] = updateRelPerm(obj)
            Sw_temp = obj.filterSw(obj.Sw);
            Krw = Sw_temp.^1.5;
            Kro = (1-Sw_temp).^1.5;
%             Krw = obj.Sw.^1.5;
%             Kro = (1-obj.Sw).^1.5;
        end
        
        %% Update derivatives of bw and bo with respect to Po.
        % [deriv_bw, deriv_bo]
        function [deriv_bw, deriv_bo] = updateDerivFVF(obj)
            % Density = exp(C*dP)
            % d(Density)/dP = C*exp(CdP)
            deriv_bw = obj.compWat*exp(obj.compWat*(obj.Po-obj.Pref));
            deriv_bo = obj.compOil*exp(obj.compOil*(obj.Po-obj.Pref));
        end
        
        %% Update derivatives of visO and visW with respect to Po.
        % [derivVisW, derivVisO]
        function [derivVisW, derivVisO] = updateDerivVis(obj)
            derivVisW = obj.compVisWat*obj.visW_std...
                *exp(obj.compVisWat*(obj.Po-obj.Pref));
            derivVisO = obj.compVisOil*obj.visO_std...
                *exp(obj.compVisOil*(obj.Po-obj.Pref));
        end
        
        %% Update derivatives of relative permeabilities with respect to Sw.
        % [derivKro, derivKrw]
        function [derivKrw, derivKro] = updateDerivRelPerm(obj)
            Sw_temp = obj.filterSw(obj.Sw);
            derivKrw = 1.5*Sw_temp.^0.5;
            derivKro = -1.5*(1-Sw_temp).^0.5;
%             derivKrw = 1.5*obj.Sw.^0.5;
%             derivKro = -1.5*(1-obj.Sw).^0.5;
        end
        
        %% Calculate derivatives of rho*g/gc with respect to Po.
        % [d(gamma_denW(iCell)/d(Po(iCell),d(gamma_den0(iCell)/d(Po(iCell)]
        function [derivDenW, derivDenO] = calcDerivDen(obj, iCell)
            derivDenW = obj.compWat*obj.gamma_denW(iCell);
            derivDenO = obj.compOil*obj.gamma_denO(iCell);
        end
        
        %% Back up properties at the previous time level.
        function [vecPPrev, PoPrev, SwPrev, bwPrev, boPrev,...
                PoPrevIter, SwPrevIter] = backUpFluid(obj)
            vecPPrev = obj.vecP;
            PoPrev = obj.Po;
            SwPrev = obj.Sw;
            boPrev = obj.bo;
            bwPrev = obj.bw;
            PoPrevIter = obj.Po;
            SwPrevIter = obj.Sw;
        end
        
        
        %% Update all fluid properties.
        function [bw, bo, Bw, Bo, gamma_denW, gamma_denO,...
                visW, visO, Krw, Kro, deriv_bw, deriv_bo,...
                derivVisW, derivVisO, derivKrw, derivKro] = updateFluid(obj)
            [bw, bo, gamma_denW, gamma_denO, Bw, Bo] = updateFVF(obj);
            [visW, visO] = updateVis(obj);
            [Krw, Kro] = updateRelPerm(obj);
            [deriv_bw, deriv_bo] = updateDerivFVF(obj);
            [derivVisW, derivVisO] = updateDerivVis(obj);
            [derivKrw, derivKro] = updateDerivRelPerm(obj);            
        end
        
        %% Filter negative Sw when calculating perm to prevent imaginary number.
        function SwNew = filterSw(obj, Sw)
            SwNew = Sw;
            for i=1:length(Sw)
                if Sw(i) < 0
                    SwNew(i) = 0;
                elseif Sw(i) > 1
                    SwNew(i) = 1;
                end
            end
        end
    end
end

