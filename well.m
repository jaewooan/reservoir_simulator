classdef well
    %Well
    %   Stores well information on each cell
    
    properties
        flowConst; % STB/day
        wellDia; % well diameter ft
        wellR; % well radius ft
        wellName; % well name
        BHPconst; %psi
        numberOfWells = 0;
        wellMap; % Map of wells(cell number->well number): 1 by nTot[prod, inj]
        wellRevMap; % Map of wells(well number->cell number): 1 by nWell[prod, inj]
        WI; % Geometric part of transmissibility[1 by nWell]
        Twell; % Well transmissibility of water and oil[2 by nWell], [water, oil]
        % Derivative of well transmissibility[4 by nWell],[dTw/dPo;dTw/dSw;dTo/dPo;dTo/dPw]
        derivTwell;
        rWB; % Approximate wellbore radius[ft][1 by nWell]
        kWB; % Wellbore absolute permeability[md][1 by nWell]
        wellType; % Prod:1, Inj:2 [1 by nWell]
        wellMode; % 1: rate, 2:BHP [1 by nWell]
        conv_alpha =0.001127; %[bbl cp/(day/psi/md/ft)]
        Pwf; % BHP[psi]
        qo; % Oil flow rate[stb/day]
        qw; % Water flow rate[stb/day]
    end
    
    methods
        function obj = well(grid)
            obj.flowConst = 0;
            obj.wellDia = 0; 
            obj.wellR = 0;
            obj.wellName = strings(1, 1);
            obj.BHPconst = 0;
            obj.wellMap = zeros(1, grid.nCell);
            obj.WI = 0;
            obj.Twell = zeros(2, 1);
            obj.derivTwell = zeros(4, 1);
            obj.wellType = 0;
            obj.rWB = 0; 
            obj.kWB = 0;
            obj.wellMode = 0;
            obj.Pwf=0;
            obj.qo=0;
            obj.qw=0;
        end  
        
        %% Calculate geometric part of well transmissibility
        function [kWB, rWB, WI] = initializeWellTransGeo(obj, grid)
            WI = obj.WI;
            kWB = 0;
            rWB = 1;
            for i = 1:obj.numberOfWells
                iCell = obj.wellRevMap(i);
                rWB(i) = 0.28*sqrt((grid.ky(iCell)./grid.kx(iCell)).^0.5*grid.dx(iCell)^2 ...
                            +(grid.kx(iCell)./grid.ky(iCell)).^0.5*grid.dy(iCell)^2)...
                             ./((grid.ky(iCell)./grid.kx(iCell)).^0.25...
                             +(grid.kx(iCell)./grid.ky(iCell)).^0.25);
                kWB(i) = sqrt(grid.kx(iCell).*grid.ky(iCell));
                WI(i) = 2*pi*obj.conv_alpha*kWB(i)*grid.dz(iCell)/log(rWB(i)/obj.wellR(i));
            end
        end
        
        %% Update well transmissibility
        function Twell = updateWellTrans(obj, fluid)
            iWellCell = obj.wellRevMap(1:obj.numberOfWells);
            phaseMobilO = 0;
            phaseMobilW = 0;
            for i = 1:length(iWellCell)
                if obj.wellType(i) == 1 % production well
                    phaseMobilW(1,i) = fluid.Krw(iWellCell(i))./...
                        (fluid.visW(iWellCell(i)).*fluid.Bw(iWellCell(i)));
                    phaseMobilO(1,i) = fluid.Kro(iWellCell(i))./...
                        (fluid.visO(iWellCell(i)).*fluid.Bo(iWellCell(i)));
                else % injection well
                    phaseMobilW(1,i) = (fluid.Krw(iWellCell(i))./fluid.visW(iWellCell(i))...
                                     +fluid.Kro(iWellCell(i))./fluid.visO(iWellCell(i)))./...
                                     fluid.Bw(iWellCell(i));
                    phaseMobilO(1,i)= fluid.Kro(iWellCell(i))./...
                        (fluid.visO(iWellCell(i)).*fluid.Bo(iWellCell(i)));
                end
            end
            Twell = [obj.WI.* phaseMobilW; obj.WI.* phaseMobilO]; %STB/day/psi            
        end
        
        %% Update derivative of well transmissiblity[dTw/dPo;dTw/dSw;dTo/dPo;dTo/dSw]
        function derivTwell = updateDerivTwell(obj, fluid)
            iWellCell = obj.wellRevMap(1:obj.numberOfWells);
            for i = 1:length(iWellCell)
                iCell = iWellCell(i);
                if obj.wellType(i) == 1 % production well
                    %dRw_i/dPo_i
                    dWatPo(1,i) = obj.WI(i).*fluid.Krw(iCell)./fluid.visW(iCell).*...
                        (fluid.deriv_bw(iCell)-fluid.bw(iCell)...
                        ./fluid.visW(iCell).*fluid.derivVisW(iCell));

                    %dRw_i/dSw_i
                    dWatSw(1,i) = obj.WI(i).*fluid.bw(iCell)...
                        ./fluid.visW(iCell).*fluid.derivKrw(iCell);

                    %dRo_i/dPo_i
                    dOilPo(1,i) = obj.WI(i).*fluid.Kro(iCell)./fluid.visO(iCell).*...
                        (fluid.deriv_bo(iCell)-fluid.bo(iCell)...
                        ./fluid.visO(iCell).*fluid.derivVisO(iCell));

                    %dRw_i/dSw_i
                    dOilSw(1,i) = obj.WI(i).*fluid.bo(iCell)./...
                        fluid.visO(iCell).*fluid.derivKro(iCell);  
                else % Injection well
                    %dRw_i/dPo_i
                    dWatPo(1,i) = obj.WI(i).*fluid.deriv_bw(iCell).*(fluid.Krw(iCell)./fluid.visW(iCell)...
                                     +fluid.Kro(iCell)./fluid.visO(iCell))+...
                                     obj.WI(i).*fluid.bw(iCell).*...
                                     (-fluid.Krw(iCell)./fluid.visW(iCell).^2*fluid.derivVisW(iCell)-...
                                     fluid.Kro(iCell)./fluid.visO(iCell).^2*fluid.derivVisO(iCell));                   
            
                    %dRw_i/dSw_i
                    dWatSw(1,i) = (fluid.derivKrw(iCell)./fluid.visW(iCell)...
                                     +fluid.derivKro(iCell)./fluid.visO(iCell)).*fluid.bw(iCell);

                    %dRo_i/dPo_i    
                    dOilPo(1,i) =0;

                    %dRw_i/dSw_i
                    dOilSw(1,i) =0;             
                end
            end         
            if length(iWellCell) == 0
                derivTwell = [0,0,0,0];
            else
                derivTwell = [dWatPo; dWatSw; dOilPo; dOilSw];
            end
        end
        
        %% Update well transmissibility and their derivatives.
        function [Twell, derivTwell] = updateWell(obj, fluid)
            Twell = obj.updateWellTrans(fluid);
            derivTwell = obj.updateDerivTwell(fluid);
        end
        
        %% Update qw, qo, and BHP
        function [qw, qo, Pwf] = updateWellCond(obj, fluid)
            for i = 1:obj.numberOfWells
                iCell = obj.wellRevMap(i);
                if obj.wellMode(i) == 1 % rate control
                    if obj.wellType(i) == 1  % production well
                        qo(i,1)=obj.flowConst(i);
                        qw(i,1)=qo(i,1)*obj.Twell(1,i)/obj.Twell(2,i);
                        Pwf(i,1)=-qo(i,1)/obj.Twell(2,i)+fluid.Po(iCell);
                    else
                        qw(i,1)= -obj.qw(1,i);%-obj.flowConst(i);
                        qo(i,1)= -obj.qo(1,i);%0;
                        Pwf(i,1)=-qw(i,1)/obj.Twell(1,i)+fluid.Po(iCell);                        
                    end
                else
                    Pwf(i,1)=obj.BHPconst(i);
                    if obj.wellType(i) == 1  % production well
                        qw(i,1)=obj.Twell(1,i)*(fluid.Po(iCell)-Pwf(i,1));
                        qo(i,1)=obj.Twell(2,i)*(fluid.Po(iCell)-Pwf(i,1));
                    else
                        qw(i,1)=obj.Twell(1,i)*(fluid.Po(iCell)-Pwf(i,1));
                        qo(i,1)=0;
                    end
                end
            end
            if obj.numberOfWells == 0
                qw = 0;
                qo = 0;
                Pwf = 0;
            end
        end
        
        %% Update qw, qo, and BHP and check if well control condition is violated
        function [isWCMOK, wellMode, Pwf, qw, qo, cumo_temp] = updateCheckWellCond(obj, fluid, dT)
            [qw, qo, Pwf] = obj.updateWellCond(fluid);
            isWCMOK = true;
            wellMode = obj.wellMode;
            cumo_temp = 0;
            for i = 1:obj.numberOfWells
                if obj.wellMode(i) == 1
                    if obj.wellType(i) == 1 && Pwf(i) < obj.BHPconst(i)
                        isWCMOK = false;
                        if abs((Pwf(i)-obj.BHPconst(i))/obj.BHPconst(i))<0.00001
                            wellMode(i) = 2;
                        end
                    elseif obj.wellType(i) == 2 && Pwf(i) > obj.BHPconst(i)
                        isWCMOK = false;
                        if abs((Pwf(i)-obj.BHPconst(i))/obj.BHPconst(i))<0.00001
                            wellMode(i) = 2;
                        end
                    end
                end
                if isWCMOK && obj.wellType(i) == 1
                    cumo_temp = cumo_temp+qo(i)*dT;
                end
            end
        end
        
        
    end
end

