classdef output
    % Saves output variables of each time
    
    properties
        outputName = cell(11,1); % Output variable names
        WOPR; % Well oil production rate [STB/day]
        WWPR; % Well water production rate [STB/day]
        WWIR; % Well water injection rate [STB/day]
        WBHPPROD; %BHP of production well [psi]
        WBHPINJ; %BHP of injection well [psi]
        FOPT; %Cumulative oil production of the field
        CFL; %CFL number
        SWAT; % Water saturation
        PORO; % Porosity
        PRESSURE; % Pressure
        numIter; % The number of iteration
        simTime;
        dT_out; % Time step(days)
        CFLmax; % Maximum CFL number
        maxR; % max norm of residual
        maxS; % max change in saturation
        maxP; % max change in pressure
    end
    
    methods
        %% Constructor
        function obj = output(grid)
            obj.WOPR = 0;
            obj.WWPR = 0;
            obj.WWIR = 0;
            obj.WBHPPROD = 0;
            obj.WBHPINJ = 0;
            obj.simTime = 0;
        end
        
        %% Add new type of variable
        function [outputName, outputValue] = addNewVar(obj, newVar)
            outputName = obj.outputName;
            outputName{1,length(outputName)+1} = newVar;
            outputValue(size(outputValue,1)+1,:,:) = zeros(1,1,size(outputValue,3));
        end
        
        %% Save values to output class
        function [WOPR, WWPR, WWIR, WBHPPROD, WBHPINJ, FOPT, SWAT, PORO, PRESSURE, CFL]...
                = updateOutput(obj, Pwf, qw, qo, cumo, well, SWAT_temp,...
                PORO_temp, PRESSURE_temp, CFL_temp)
            iInj = 0; iProd = 0;
            for i = 1:well.numberOfWells
                if well.wellType(i) == 1 % production
                    iProd = iProd + 1;
                    WOPR(iProd,:)=qo(i,:);
                    WWPR(iProd,:)=qw(i,:);
                    WBHPPROD(iProd,:) = Pwf(i,:);
                else % injection
                    iInj = iInj + 1;
                    WWIR(iInj,:) = qw(i,:);   
                    WBHPINJ(iInj,:) = Pwf(i,:);   
                end
            end
            if well.numberOfWells == 0
                WOPR = 0;
                WWPR = 0;
                WBHPPROD = 0;
                WWIR = 0; 
                WBHPINJ = 0;
            end
            FOPT = cumo;
            PORO = PORO_temp;
            SWAT = SWAT_temp;
            PRESSURE = PRESSURE_temp;
            CFL = CFL_temp;
%             for row = 1:20
%                 for col = 1:20
%                     PORO(:,row,col) = PORO_temp(:,col+20*(row-1));
%                     SWAT(:,row,col) = SWAT_temp(:,col+20*(row-1));
%                     PRESSURE(:,row,col) = PRESSURE_temp(:,col+20*(row-1));
%                     CFL(:,row,col) = CFL_temp(:,col+20*(row-1));
%                 end
%             end
        end
    end
end

