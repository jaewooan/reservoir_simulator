classdef equation
    % Class that stores variables and functions related to Newton method.
    % Generate residual and Jacobian.
    % Solve the problem with Newton method.
    %%
    properties
        %% Information for sparse matrix
        iDim; % Dimension in x direction
        jDim; % Dimension in y direction
        kDim; % Dimension in z direction
        nDim; % Total dimension of this problem
        matrixJ; % Jacobian matrix
        arrayR; % Residual matrix.
        iRowNonZeroJMatrix; % Row index of Jacobian nonzero element
        iColNonZeroJMatrix; % Column index of Jacobian nonzero element
        % Cell index corresponding to each row of Jacobian matrix
        iCellRowNonZeroJMatrix; 
        % Cell index corresponding to cach column of Jacobian matrix
        iCellColNonZeroJMatrix; 
        % Index of face related to each off-diagonal matrix component
        iFaceOfMatrix; 
        % Index of face related to each diagonal matrix component
        iFaceDiag;         
        %% Timestepping and convergence controls(default setting in Eclipse).
        epsR = 1e-3;% Maximum allowable residual in convergence
        epsS = 1e-3;% Maximum allowable saturation change in convergence
        epsP = 1e-3;% Maximum allowable pressure change in convergence        
        %% Iteration-related variables        
        delT = 0.001; % Simulation time step(days)
        curT = 0; % Current simulation time(days)
        totT = 0; % Total simulation time(days)
        iterT = 0; % The number of iteration over time loop.
        iterK = 0; % The number of iteration over Newton iteration
        maxIterK = 10; % The maximum number of iteration over Newton-Rhapson
    end
    methods
        
        function obj = equation(grid, connection)
            nTot = grid.nCell;
            nTotMatrix = grid.nMatrix;
            nComp = 1; % The number of matrix component
            obj.arrayR = zeros(1, 2*nTot);
            obj.matrixJ = zeros(2*nTot, 2*nTot);
            obj.iFaceDiag = cell(1,nTot);
            % Initialize the location of Jacobian component which is used
            % to make sparse matrix; iRowNonZeroJMatrix and
            % iColNonZeroJMatrix.
            for iCell=1:nTot    
                iCellNeighbors = connection.iNeighbor{iCell};
                vCellNonzero = [iCellNeighbors,iCell];
                for j=1:length(vCellNonzero)
                    jCell = vCellNonzero(j);
                    % Remember the location of nonzero matrix element
                    obj.iRowNonZeroJMatrix(nComp:nComp+3) = [2*iCell-1,2*iCell,2*iCell-1,2*iCell];
                    obj.iColNonZeroJMatrix(nComp:nComp+3) = [2*jCell-1,2*jCell-1,2*jCell,2*jCell];
                    % Remember the cell corresponding to row and column of matrix element
                    obj.iCellRowNonZeroJMatrix(nComp:nComp+3) = [iCell,iCell,iCell,iCell];
                    obj.iCellColNonZeroJMatrix(nComp:nComp+3) = [jCell,jCell,jCell,jCell];
                    if iCell~=jCell
                        % Remember the face index related to the off-diagonal matrix
                        % element
                        obj.iFaceOfMatrix(nComp:nComp+3) = connection.iFace(iCell,jCell)*ones(1,4);
                    end
                    % Move to the next 2 by 2 block matrix element.
                    nComp = nComp+4;
                end
                % Remember the face index related to the diagonal matrix element
                obj.iFaceDiag(iCell) = {connection.iFace(iCell,iCellNeighbors)}; 
            end 
        end  
        
        %% Generate class
%         function obj = equation(grid, connection)
%             nTot = grid.nCell;
%             nTotMatrix = grid.nMatrix;
%             nComp = 1; % The number of matrix component
%             obj.arrayR = zeros(1, 2*nTot);
%             obj.matrixJ = zeros(2*nTot, 2*nTot);
%             obj.iFaceDiag = cell(1,nTot);
%             % The dimension-related number for checking the neighbors of
%             % each cell.
%             obj.iDim = (grid.nx>1); obj.jDim = (grid.ny>1);
%             obj.kDim = (grid.nz>1); obj.nDim = obj.iDim+obj.jDim+obj.kDim;
%             % Initialize the location of Jacobian component which is used
%             % to make sparse matrix; iRowNonZeroJMatrix and
%             % iColNonZeroJMatrix.
%             for i=1:nTotMatrix
%                 for j=1:nTotMatrix
%                             % Neighbor in -z direction
%                     if (obj.kDim == 1 && j == i-grid.nx*grid.ny && j>=1) ||... 
%                          ...% Neighbor in -y direction
%                         (obj.jDim == 1 && j == i-grid.nx && ...
%                             floor((i-1)/(grid.nx*grid.ny))==floor((j-1)/(grid.nx*grid.ny))) ||...
%                          ...% Neighbor in -x direction
%                         (obj.iDim == 1 && j == i-1 && floor((i-1)/grid.nx)==floor((j-1)/grid.nx)) ||... 
%                         (i==j) ||... % Itself
%                          ...% Neighbor in +x direction
%                         (obj.iDim == 1 && j == i+1 && floor((i-1)/grid.nx)==floor((j-1)/grid.nx)) ||... 
%                          ...% Neighbor in +y direction
%                         (obj.jDim == 1 && j == i+grid.nx && ...
%                             floor((i-1)/(grid.nx*grid.ny))==floor((j-1)/(grid.nx*grid.ny))) ||...
%                          ...% Neighbor in +z direction
%                         (obj.kDim == 1 && j == i+grid.nx*grid.ny && j<=grid.nx*grid.ny*grid.nz)
%                         % Remember the location of nonzero matrix element
%                         obj.iRowNonZeroJMatrix(nComp:nComp+3) = [2*i-1,2*i,2*i-1,2*i];
%                         obj.iColNonZeroJMatrix(nComp:nComp+3) = [2*j-1,2*j-1,2*j,2*j];
%                         % Remember the cell corresponding to row and column of matrix element
%                         obj.iCellRowNonZeroJMatrix(nComp:nComp+3) = [i,i,i,i];
%                         obj.iCellColNonZeroJMatrix(nComp:nComp+3) = [j,j,j,j];
%                         if i~=j
%                             % Remember the face index related to the off-diagonal matrix
%                             % element
%                             obj.iFaceOfMatrix(nComp:nComp+3) = connection.iFace(i,j)*ones(1,4);
%                         end
%                         % Move to the next 2 by 2 block matrix element.
%                         nComp = nComp+4;
%                     end                         
%                 end
%                 % Remember the face index related to the diagonal matrix element
%                 obj.iFaceDiag(i) = {connection.iFace(i,connection.iNeighbor{i})}; 
%             end 
%         end        
        %% Make residual array
        function [arrayR, q_out] = makeResidual(obj, grid, fluid, connection, well)
            nTot =grid.nCell; %grid.nx*grid.ny*grid.nz; 
            % For each residual i,
            for i = 1:nTot
                % Read the faces and neighboring cells related to residual i
                iFace = connection.iFace(i,connection.iNeighbor{i});
                iNeig = connection.iNeighbor{i};
                flowPart_water=0;
                flowPart_oil=0;
                flowPart_water_out=0;
                flowPart_oil_out=0;
                if i==442 || sum(ismember(iNeig, 442))
                    %disp('dd')
                end
                % For each face,
                for j=1:length(iFace)
                    % Sum the flow terms in oil and water residuals
                    try
                    flowPart_water_temp = connection.transTot(1,iFace(j))...
                        .*(fluid.Po(iNeig(j))-fluid.Po(i)...
                        -connection.iGravity(iFace(j)).*connection.avgDensity(1,iFace(j))...
                        .*(grid.zCell(iNeig(j))-grid.zCell(i)));
                    catch
                        disp('dd')
                    end
                    flowPart_oil_temp = connection.transTot(2,iFace(j))...
                        .*(fluid.Po(iNeig(j))-fluid.Po(i)...
                        -connection.iGravity(iFace(j)).*connection.avgDensity(2,iFace(j))...
                        .*(grid.zCell(iNeig(j))-grid.zCell(i)));
                    flowPart_water = flowPart_water+flowPart_water_temp;
                    flowPart_oil = flowPart_oil+flowPart_oil_temp;
                    flowPart_water_out = flowPart_water_out + connection.iUpwind(1,iFace(j))*flowPart_water_temp*fluid.Bw(i); %RB/day
                    flowPart_oil_out = flowPart_oil_out + connection.iUpwind(2,iFace(j))*flowPart_oil_temp*fluid.Bo(i); %RB/day
                end
                % Calculate the accumulation term in oil and water
                % residuals.
                accumPart_water=grid.volCell(i)./obj.delT/5.615...
                          .*(grid.poro(i)*fluid.Sw(i)*fluid.bw(i)...
                             -grid.poroPrev(i)*fluid.SwPrev(i)*fluid.bwPrev(i));
                accumPart_oil=grid.volCell(i)/obj.delT/5.615...
                          .*(grid.poro(i)*(1-fluid.Sw(i))*fluid.bo(i)...
                             -grid.poroPrev(i)*(1-fluid.SwPrev(i))*fluid.boPrev(i));                         
                % Calculate well part
                wellPart_water=0;
                wellPart_oil=0;
                iWell = well.wellMap(1,i);
                if iWell > 0 % If well exists,                
                    if well.wellMode(1,iWell) == 1 %Rate control
                        if well.wellType(iWell) == 1 % production well
                            wellPart_water = well.flowConst(iWell)*well.Twell(1,iWell)/well.Twell(2,iWell);
                            wellPart_oil = well.flowConst(iWell);
                        else % injection well
                            wellPart_water = -well.flowConst(iWell);%--well.qw(1,iWell)
                            wellPart_oil = 0;%-well.qo(1,iWell);%0;
                        end
                    else %BHP control
                        if well.wellType(iWell) == 1 % production well
                            wellPart_water = well.Twell(1,iWell)*(fluid.Po(i)-well.BHPconst(iWell));
                            wellPart_oil = well.Twell(2,iWell)*(fluid.Po(i)-well.BHPconst(iWell));
                        else
                            wellPart_water = well.Twell(1,iWell)*(fluid.Po(i)-well.BHPconst(iWell));
                            wellPart_oil = 0;
                        end
                    end                                   
                end                
                flowPart_oil_bc = connection.bcTerm(i);
                % Calculate residual: flow term - accumulation term
                arrayR(2*i-1) = flowPart_water-accumPart_water-wellPart_water;
                arrayR(2*i) =flowPart_oil-accumPart_oil-wellPart_oil+flowPart_oil_bc;
                q_out(i) = abs(flowPart_water_out+flowPart_oil_out);
            end
            arrayR = sparse(arrayR');
        end        
        %% Make Jacobian matrix
        % Structure of variables.
        % transTot(2 by nFace): [Tw, To]
        
        % avgDensity(2 by nFace):
        % [(gamma_denW(iCell)+gamma_denW(iCellNext)/2,
        %  (gamma_denO(iCell)+gamma_denO(iCellNext)/2]
        
        % derivTransTot(8 by nFace):   
        %[dTw(iFace)/dPo(iCell), dTw(iFace)/dPo(iCellNeighbor),
        % dTw(iFace)/dSw(iCell), dTw(iFace)/dSw(iCellNeighbor),
        % dTo(iFace)/dPo(iCell), dTo(iFace)/dPo(iCellNeighbor),
        % dTo(iFace)/dSw(iCell), dTo(iFace)/dSw(iCellNeighbor)]
        
        % derivAvgDensity(4 by nFace):
        %[dDenW(iFace)/dPo(iCell), dDenW(iFace)/dPo(iCellNext),
        % dDenO(iFace)/dPo(iCell), dDenO(iFace)/dPo(iCellNext)]
        function matrixJ = makeJacobian(obj, grid, fluid, connection, well)
            compMatrixJ = zeros(1,length(obj.iRowNonZeroJMatrix));
            % Calculate [dRw/dPo, dRo/dPo, dRw/dSw, dRo/dSw] of each block 
            % matrix at each iteration
            for nComp = 1:4:length(obj.iCellRowNonZeroJMatrix)
                % i: the cell index corresponding to the ith oil and
                % water residuals(ex: dRo_i/dP_j)
                i = obj.iCellRowNonZeroJMatrix(nComp);             
                % Off-diagonal block
                if  i~= obj.iCellColNonZeroJMatrix(nComp)
                    % Read the cell index of Po and Sw in the denominator
                    % of partial derivative of the block diagonal(ex: dRo_i/dP_j)
                    j = obj.iCellColNonZeroJMatrix(nComp);
                    iFace = obj.iFaceOfMatrix(nComp);
                    delZ = grid.zCell(j)-grid.zCell(i);
                    delP = fluid.Po(j)-fluid.Po(i);                    
                    %dRw_i/dPo_j & dRo_i/dPo_j
                    compMatrixJ(nComp:nComp+1)=...
                        connection.transTot(1:2, iFace).*(...
                            1-connection.iGravity(iFace).*...
                            connection.derivAvgDensity([2,4],iFace).*delZ)...
                       +connection.derivTransTot([2,6],iFace).*(...
                             delP-connection.iGravity(iFace).*...
                             connection.avgDensity([1,2],iFace).*delZ);
                    %dRw_i/dSw_j, dRo_i/dSw_j
                    compMatrixJ(nComp+2:nComp+3) =...
                        connection.derivTransTot([4,8],iFace).*...
                            (delP-connection.iGravity(iFace).*...
                                connection.avgDensity([1,2],iFace).*delZ);
                else % Diagonal block
                    iFace = obj.iFaceDiag{i};
                    j = connection.iNeighbor{i};
                    delZ = grid.zCell(j)-grid.zCell(i);
                    delP = fluid.Po(j)-fluid.Po(i);                    
                    %Calculate Well parts                    
                    dqwdpo=0; dqwdSw=0; dqodpo=0; dqodSw=0;
                    iWell = well.wellMap(1,i);
                    if iWell > 0 % If well exists,   
                        %[dTw/dPo;dTw/dSw;dTo/dPo;dTo/dSw]
                        if well.wellMode(1,iWell) == 1 %Rate control
                            if well.wellType(iWell) == 1% production
                                dqwdpo = well.flowConst(iWell)*(...
                                    well.derivTwell(1,iWell)/well.Twell(2,iWell)...
                                    -well.Twell(1,iWell)/well.Twell(2,iWell)^2*...
                                    well.derivTwell(3,iWell));
                                dqwdSw = well.flowConst(iWell)*(...
                                    well.derivTwell(2,iWell)/well.Twell(2,iWell)...
                                    -well.Twell(1,iWell)/well.Twell(2,iWell)^2*...
                                    well.derivTwell(4,iWell));
                            end
                        else %BHP control
                            if well.wellType(iWell) == 1 % production
                                dqwdpo=well.derivTwell(1,iWell)*...
                                    (fluid.Po(i)-well.BHPconst(iWell))...
                                    +well.Twell(1,iWell);
                                dqwdSw=well.derivTwell(2,iWell)*...
                                    (fluid.Po(i)-well.BHPconst(iWell));
                                dqodpo=well.derivTwell(3,iWell)*...
                                    (fluid.Po(i)-well.BHPconst(iWell))...
                                    +well.Twell(2,iWell);
                                dqodSw=well.derivTwell(4,iWell)*...
                                    (fluid.Po(i)-well.BHPconst(iWell));
                            else
                                dqwdpo=well.derivTwell(1,iWell)*...
                                    (fluid.Po(i)-well.BHPconst(iWell))...
                                    +well.Twell(1,iWell);
                                dqwdSw=well.derivTwell(2,iWell)*...
                                    (fluid.Po(i)-well.BHPconst(iWell));
                            end
                        end                                   
                    end                    
                    %dRw_i/dPo_j
                    compMatrixJ(nComp)=...
                        sum(connection.derivTransTot(1,iFace).*...
                            (delP-connection.iGravity(iFace).*...
                             connection.avgDensity(1,iFace).*delZ))...
                       -sum(connection.transTot(1,iFace).*...
                            (1+connection.iGravity(iFace).*...
                              connection.derivAvgDensity(1,iFace).*delZ))...
                       -grid.volCell(i)*fluid.Sw(i)/5.615/obj.delT.*...
                           (grid.poro(i).*fluid.deriv_bw(i)+...
                            fluid.bw(i).*grid.derivPoro(i))-dqwdpo;                        
                    %dRo_i/dPo_j
                    compMatrixJ(nComp+1)=...
                        sum(connection.derivTransTot(5,iFace).*...
                            (delP-connection.iGravity(iFace).*...
                             connection.avgDensity(2,iFace).*delZ))...
                       -sum(connection.transTot(2,iFace).*...
                            (1+connection.iGravity(iFace).*...
                              connection.derivAvgDensity(3,iFace).*delZ))...
                       -grid.volCell(i)*(1-fluid.Sw(i))/5.615/obj.delT.*...
                           (grid.poro(i).*fluid.deriv_bo(i)+...
                            fluid.bo(i).*grid.derivPoro(i))-dqodpo; 
                    compMatrixJ(nComp+1) = compMatrixJ(nComp+1) + connection.bcDerivTerm(i);% add bc
                    %dRw_i/dSw_j
                    compMatrixJ(nComp+2) =...
                        sum(connection.derivTransTot(3,iFace).*...
                            (delP-connection.iGravity(iFace).*...
                             connection.avgDensity(1,iFace).*delZ))...
                       -grid.volCell(i)/5.615/obj.delT.*...
                           (grid.poro(i).*fluid.bw(i))-dqwdSw;                       
                    %dRo_i/dSw_j
                    compMatrixJ(nComp+3) =...
                        sum(connection.derivTransTot(7,iFace).*...
                            (delP-connection.iGravity(iFace).*...
                             connection.avgDensity(2,iFace).*delZ))...
                       +grid.volCell(i)/5.615/obj.delT.*...
                           (grid.poro(i).*fluid.bo(i))-dqodSw;
                end
                % Move to the next 2 by 2 block matrix
                nComp = nComp + 4;
            end
            matrixJ = sparse(obj.iRowNonZeroJMatrix, ...
                obj.iColNonZeroJMatrix, compMatrixJ);
        end
        
        %% Update time step
        function dT = updateTimeStep(obj, fluid, isWCMOK)
            dXs_max = 0.1*max(abs(fluid.Sw));
            dXp_max = 0.01*max(abs(fluid.Po));            
            w = 0.5;            
            if ~isWCMOK 
                dT = obj.delT; 
            else
                dXs = fluid.Sw-fluid.SwPrev;
                dXp = fluid.Po-fluid.PoPrev;
                dT = obj.delT*(1+w)*min([dXs_max./(abs(dXs)+w*dXs_max),...
                    dXp_max./(abs(dXp)+w*dXp_max)]);
            end            
            if dT + obj.curT > obj.totT
                dT = obj.totT - obj.curT;
            end            
            if dT >= 100 % Maximum time limit
                dT = 100;
            end
        end        
        %% Solve the problem using Newton-Rhapson method
        function [isConv, vecP, Po, Sw, bo, bw, poro, well, k, q_out,...
                maxR, maxS, maxP] = solveNRProblem(obj, fluid, grid, connection, well)
            k = 0;
            isConv = false;
            R = zeros(1, 2);
            while ~isConv
                % Is k max iteration?
                if k >= obj.maxIterK
                    isConv = false;
                    break;
                end
                k = k + 1;
                
                % Update all properties
                [fluid.bw, fluid.bo, fluid.Bw, fluid.Bo, fluid.gamma_denW, fluid.gamma_denO,...
                 fluid.visW, fluid.visO, fluid.Krw, fluid.Kro, fluid.deriv_bw, fluid.deriv_bo,...
                 fluid.derivVisW, fluid.derivVisO, fluid.derivKrw, fluid.derivKro]...
                 = fluid.updateFluid();
                [grid.poro, grid.derivPoro] = grid.updateGrid(fluid);
                [connection.iUpwind, connection.transTot, connection.avgDensity,...
                 connection.derivTransTot, connection.derivAvgDensity, ...
                 connection.bcTerm, connection.bcDerivTerm]...
                 = connection.updateConnnection(fluid, grid);  
                [well.Twell, well.derivTwell] = well.updateWell(fluid);
                % Make Jacobian and residual.
                [R, q_out] = obj.makeResidual(grid, fluid, connection, well);
                J = obj.makeJacobian(grid, fluid, connection, well);
                % Update pressure and saturation.
                dX = -(J\R)';
                realIndex = abs(imag(dX)) < 1e-10;
                dX(realIndex) = real(dX(realIndex));
                fluid.PoPrevIter = fluid.Po; fluid.SwPrevIter = fluid.Sw;
                fluid.vecP = fluid.vecP + dX;
                fluid.Po = fluid.vecP(1:2:end); fluid.Sw = fluid.vecP(2:2:end);
                [isConv, maxR, maxS, maxP] = obj.checkConvergence(fluid, grid, R, k);
            end            
            % Result   vecP, Po, Sw, bo, bw, poro
            vecP = fluid.vecP;
            Po = fluid.Po;
            Sw = fluid.Sw;
            [bw, bo] = fluid.updatebwbo(Po);
            poro = grid.updatePoro(Po); 
        end        
        %% Check convergence
        function [isConv, maxR, maxS, maxP] = checkConvergence(obj, fluid, grid, arrayR, k)
            % Infinity norm of the residual
            maxRw = norm(abs(5.615*fluid.Bw.*arrayR(1:2:end)*obj.delT./...
                            (grid.volCell.*grid.poro)), Inf);
            maxRo = norm(abs(5.615*fluid.Bo.*arrayR(2:2:end)*obj.delT./...
                            (grid.volCell.*grid.poro)), Inf);
            maxR = max(maxRw, maxRo);
            % Max change in saturation
            maxS = norm(abs(fluid.Sw - fluid.SwPrevIter), Inf);
            % Max change in pressure
            maxP = norm(abs((fluid.Po - fluid.PoPrevIter)/mean(fluid.Po)), Inf);            
            if maxR < obj.epsR && maxS < obj.epsS && maxP < obj.epsP && k > 0
                isConv = true;
            else
                isConv = false;
            end
        end
    end
end

