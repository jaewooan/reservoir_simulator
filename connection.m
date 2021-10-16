classdef connection
    %CELLCONNECTION: Store faces connected to a cell.
    properties
        iCell; % Cell index of each face [1 by nFace]
        iCellNext; % Neighbor cell index to iCell in each face [1 by nFace]
        iNeighbor; % Neighbors to each cell.
        iFace; % Face index of each face from cell i to j
        iGravity; % Indicator there is gravity effect across the face
        iUpwind; % Indicator if potential of iCell is larger than that of iCellNext
        nFace; % The number of faces.
        areaFace; % Area[ft2] of each cell [1 by nFace]
        distFace; % Cell width[ft] of two cells [2 by nFace]
        permFace; % Absolute permeability[md] of two cells [2 by nFace]
        conv_alpha =0.001127; %[bbl cp/(day/psi/md/ft)]
        % Geometric part of transmissibility of each face [1 by nFace]
        transGeo;
        % Total transmissibility of each face. [2 by nFace] [Tw; To]
        transTot;
        % Average density of densities of cell i and j in each face
        % [2 by nFace]
        % [(gamma_denW(iCell)+gamma_denW(iCellNext)/2;
        %  (gamma_denO(iCell)+gamma_denO(iCellNext)/2]
        avgDensity;        
        % Derivatives of transTot with respect to Po ans Sw
        % [8 by nFace]
        %[dTw(iFace)/dPo(iCell); dTw(iFace)/dPo(iCellNeighbor);
        % dTw(iFace)/dSw(iCell); dTw(iFace)/dSw(iCellNeighbor);
        % dTo(iFace)/dPo(iCell); dTo(iFace)/dPo(iCellNeighbor);
        % dTo(iFace)/dSw(iCell); dTo(iFace)/dSw(iCellNeighbor)]
        derivTransTot;   
        % Derivatives of average densities with respect to Po
        % derivAvgDensity: 4 by nFace matrix
        %[dDenW(iFace)/dPo(iCell); dDenW(iFace)/dPo(iCellNext);
        % dDenO(iFace)/dPo(iCell); dDenO(iFace)/dPo(iCellNext)]
        derivAvgDensity;
        % The type of connections
        % 0: matrix - matrix
        % 1: matrix - fracture
        % 2: fracture - fracture
        iConnectType;
        % PL, PR
        PL;
        PR;
        bcTerm;
        bcDerivTerm;
    end    
    methods
        %% Class initialization
        function obj = connection(grid)
            % Initialize conneciton
            iFaceTemp = 1;
            obj.iNeighbor = cell(1,grid.nx*grid.ny*grid.nz);
            for k = 1:grid.nz
                for j = 1:grid.ny
                    for i = 1:grid.nx
                        iCell_1D = i + (j-1)*grid.nx + (k-1)*grid.nx*grid.ny;                        
                        if i < grid.nx % a right neighbor cell (i+1, j, k)
                            iCellNext_1D = iCell_1D + 1;
                            obj.iFace(iCell_1D, iCellNext_1D) = iFaceTemp;
                            obj.iCell(iFaceTemp) = iCell_1D;
                            obj.iCellNext(iFaceTemp) = iCellNext_1D;
                            obj.iGravity(iFaceTemp) = 0;
                            obj.areaFace(iFaceTemp) = grid.dy(1,iCell_1D)*grid.dz(1,iCell_1D);
                            obj.distFace(:, iFaceTemp) = [grid.dx(1,iCell_1D), grid.dx(1,iCell_1D)];
                            obj.permFace(:, iFaceTemp) = [grid.kx(iCell_1D), grid.kx(iCellNext_1D)];
                            obj.iNeighbor(iCell_1D) = {[obj.iNeighbor{iCell_1D},iCellNext_1D]};
                            obj.iConnectType(1,iFaceTemp) = 0;
                            iFaceTemp = iFaceTemp + 1;
                        end
                        if i > 1  % a left neighbor cell (i-1, j, k)
                            iCellNext_1D = iCell_1D - 1;
                            obj.iFace(iCell_1D, iCellNext_1D) = iFaceTemp;
                            obj.iCell(iFaceTemp) = iCell_1D;
                            obj.iCellNext(iFaceTemp) = iCellNext_1D;
                            obj.iGravity(iFaceTemp) = 0;
                            obj.areaFace(iFaceTemp) = grid.dy(1,iCell_1D)*grid.dz(1,iCell_1D);
                            obj.distFace(:, iFaceTemp) = [grid.dx(1,iCell_1D), grid.dx(1,iCell_1D)];
                            obj.permFace(:, iFaceTemp) = [grid.kx(iCell_1D), grid.kx(iCellNext_1D)];
                            obj.iNeighbor(iCell_1D) = {[obj.iNeighbor{iCell_1D},iCellNext_1D]};
                            obj.iConnectType(1,iFaceTemp) = 0;
                            iFaceTemp = iFaceTemp + 1;
                        end
                        if j < grid.ny  % a front neighbor cell (i, j+1, k)                            
                            iCellNext_1D = iCell_1D + grid.nx;
                            obj.iFace(iCell_1D, iCellNext_1D) = iFaceTemp;
                            obj.iCell(iFaceTemp) = iCell_1D;
                            obj.iCellNext(iFaceTemp) = iCellNext_1D;
                            obj.iGravity(iFaceTemp) = 0;
                            obj.areaFace(iFaceTemp) = grid.dx(1,iCell_1D)*grid.dz(1,iCell_1D);
                            obj.distFace(:, iFaceTemp) = [grid.dy(1,iCell_1D), grid.dy(1,iCell_1D)];
                            obj.permFace(:, iFaceTemp) = [grid.ky(iCell_1D), grid.ky(iCellNext_1D)];
                            obj.iNeighbor(iCell_1D) = {[obj.iNeighbor{iCell_1D},iCellNext_1D]};
                            obj.iConnectType(1,iFaceTemp) = 0;
                            iFaceTemp = iFaceTemp + 1;
                        end
                        if j > 1  % a back neighbor cell (i, j-1, k)                            
                            iCellNext_1D = iCell_1D - grid.nx;
                            obj.iFace(iCell_1D, iCellNext_1D) = iFaceTemp;
                            obj.iCell(iFaceTemp) = iCell_1D;
                            obj.iCellNext(iFaceTemp) = iCellNext_1D;
                            obj.iGravity(iFaceTemp) = 0;
                            obj.areaFace(iFaceTemp) = grid.dx(1,iCell_1D)*grid.dz(1,iCell_1D);
                            obj.distFace(:, iFaceTemp) = [grid.dy(1,iCell_1D), grid.dy(1,iCell_1D)];
                            obj.permFace(:, iFaceTemp) = [grid.ky(iCell_1D), grid.ky(iCellNext_1D)];
                            obj.iNeighbor(iCell_1D) = {[obj.iNeighbor{iCell_1D},iCellNext_1D]};
                            obj.iConnectType(1,iFaceTemp) = 0;
                            iFaceTemp = iFaceTemp + 1;
                        end
                        if k < grid.nz  % a upper neighbor cell (i, j, k+1)                                                        
                            iCellNext_1D = iCell_1D + grid.nx*grid.ny;
                            obj.iFace(iCell_1D, iCellNext_1D) = iFaceTemp;
                            obj.iCell(iFaceTemp) = iCell_1D;
                            obj.iCellNext(iFaceTemp) = iCellNext_1D;
                            obj.iGravity(iFaceTemp) = 1;
                            obj.areaFace(iFaceTemp) = grid.dx(1,iCell_1D)*grid.dy(1,iCell_1D);
                            obj.distFace(:, iFaceTemp) = [grid.dz(1,iCell_1D), grid.dz(1,iCell_1D)];
                            obj.permFace(:, iFaceTemp) = [grid.kz(iCell_1D), grid.kz(iCellNext_1D)];
                            obj.iNeighbor(iCell_1D) = {[obj.iNeighbor{iCell_1D},iCellNext_1D]};
                            obj.iConnectType(1,iFaceTemp) = 0;
                            iFaceTemp = iFaceTemp + 1;
                        end
                        if k > 1  % a lower neighbor cell (i, j, k-1)                                                       
                            iCellNext_1D = iCell_1D - grid.nx*grid.ny;
                            obj.iFace(iCell_1D, iCellNext_1D) = iFaceTemp;
                            obj.iCell(iFaceTemp) = iCell_1D;
                            obj.iCellNext(iFaceTemp) = iCellNext_1D;
                            obj.iGravity(iFaceTemp) = 1;
                            obj.areaFace(iFaceTemp) = grid.dx(1,iCell_1D)*grid.dy(1,iCell_1D);
                            obj.distFace(:, iFaceTemp) = [grid.dz(1,iCell_1D), grid.dz(1,iCell_1D)];
                            obj.permFace(:, iFaceTemp) = [grid.kz(iCell_1D), grid.kz(iCellNext_1D)];
                            obj.iNeighbor(iCell_1D) = {[obj.iNeighbor{iCell_1D},iCellNext_1D]};
                            obj.iConnectType(1,iFaceTemp) = 0;
                            iFaceTemp = iFaceTemp + 1;
                        end
                    end
                end
            end
            obj.nFace = iFaceTemp-1;            
            % Build vacant matrices.
            obj.transGeo = zeros(1, obj.nFace);
            obj.transTot = zeros(2, obj.nFace);
            obj.avgDensity = zeros(2, obj.nFace);
            obj.derivTransTot = zeros(8, obj.nFace);
            obj.derivAvgDensity = zeros(4, obj.nFace);
        end    
        
        %% Shift arrays
        function obj = shiftArray(obj, ID_MCell_1, ID_MCell_2) 
            nxSize = size(obj.iFace,1);
            nySize = size(obj.iFace,2);                      
            iPrev = 1:ID_MCell_1;
            iNext = ID_MCell_2:length(obj.iNeighbor);
            
            % Shift arrays          
            obj.iFace((ID_MCell_2+1):(nxSize+1), :) = obj.iFace(ID_MCell_2:nxSize, :);
            obj.iFace(:, (ID_MCell_2+1):(nySize+1)) = obj.iFace(:, ID_MCell_2:nySize);  
            obj.iNeighbor = [obj.iNeighbor(1,iPrev), {[]}, obj.iNeighbor(1,iNext)];
            
            % Make vacant spaces zero
            obj.iFace(ID_MCell_2, :) = 0;
            obj.iFace(:, ID_MCell_2) = 0;
            
            % Increase by one
            for i=1:length(obj.iCell)
                if obj.iCell(i) >= ID_MCell_2
                    obj.iCell(i) = obj.iCell(i) + 1;
                end
                
                if obj.iCellNext(i) >= ID_MCell_2
                    obj.iCellNext(i) = obj.iCellNext(i) + 1;                    
                end
            end
            
            for i=1:length(obj.iNeighbor)
                vNeighbor = obj.iNeighbor{i};
                iIncrease = (vNeighbor>=ID_MCell_2);
                vNeighbor(iIncrease) = vNeighbor(iIncrease) + 1;
                obj.iNeighbor(i) = {vNeighbor};
            end
        end
        
        %% Copy connection of M1 to M2.
        function obj = copyConnection(obj, grid, ID_MCell_1, ID_MCell_2)        
            vNeighbor = obj.iNeighbor{ID_MCell_1};
            for i = 1:length(vNeighbor)                
                % Split matrix -> neighbor
                obj.nFace = obj.nFace + 1;
                iFaceTemp = obj.nFace;
                iCell_1D = ID_MCell_2;
                iCellNext_1D = vNeighbor(i);
                iFaceRep = obj.iFace(ID_MCell_1, iCellNext_1D);
                obj.iFace(iCell_1D, iCellNext_1D) = iFaceTemp;
                obj.iCell(iFaceTemp) = iCell_1D;
                obj.iCellNext(iFaceTemp) = iCellNext_1D;
                obj.iGravity(iFaceTemp) = 0;
                obj.areaFace(iFaceTemp) = obj.areaFace(iFaceRep);
                obj.distFace(:, iFaceTemp) = obj.distFace(:, iFaceRep);
                obj.permFace(:, iFaceTemp) = [grid.kx(iCell_1D), grid.kx(iCellNext_1D)];
                if length(obj.iNeighbor{iCell_1D}) > 0
                    obj.iNeighbor(iCell_1D) = {[obj.iNeighbor{iCell_1D},iCellNext_1D]};
                else
                    obj.iNeighbor(iCell_1D) = {[iCellNext_1D]};                    
                end
                obj.iConnectType(1,iFaceTemp) = 0;                
                                
                % neighbor -> Split matrix 
                obj.nFace = obj.nFace + 1;
                iFaceTemp = obj.nFace;
                iCell_1D = vNeighbor(i);
                iCellNext_1D = ID_MCell_2;
                iFaceRep = obj.iFace(iCell_1D, ID_MCell_1);
                obj.iFace(iCell_1D, iCellNext_1D) = iFaceTemp;
                obj.iCell(iFaceTemp) = iCell_1D;
                obj.iCellNext(iFaceTemp) = iCellNext_1D;
                obj.iGravity(iFaceTemp) = 0;
                obj.areaFace(iFaceTemp) = obj.areaFace(iFaceRep);
                obj.distFace(:, iFaceTemp) = obj.distFace(:, iFaceRep);
                obj.permFace(:, iFaceTemp) = [grid.kx(iCell_1D), grid.kx(iCellNext_1D)];
                if length(obj.iNeighbor{iCell_1D}) > 0
                    obj.iNeighbor(iCell_1D) = {[obj.iNeighbor{iCell_1D},iCellNext_1D]};
                else
                    obj.iNeighbor(iCell_1D) = {[iCellNext_1D]};                    
                end
                obj.iConnectType(1,iFaceTemp) = 0;
            end
        end
        
        
         %% Delete connections
        function obj = deleteConnection(obj, iCell1, iCell2)            
            for i=1:2
                if i==1
                    iCell_1D = iCell1;
                    iCellNext_1D = iCell2;
                else
                    iCell_1D = iCell2;
                    iCellNext_1D = iCell1;
                end
                obj.nFace = obj.nFace - 1;
                iFaceTemp = obj.iFace(iCell_1D, iCellNext_1D);
                obj.iFace(iCell_1D, iCellNext_1D) = 0;
                obj.iCell(iFaceTemp) = [];
                obj.iCellNext(iFaceTemp) = [];
                obj.iGravity(iFaceTemp) = [];
                obj.areaFace(iFaceTemp) = [];
                obj.distFace(:, iFaceTemp) = [];
                obj.permFace(:, iFaceTemp) = [];
                vNeighbor = obj.iNeighbor{iCell_1D};
                for j=1:length(vNeighbor)
                    if vNeighbor(j) == iCellNext_1D
                        vNeighbor(j) = [];
                        break;
                    end
                end            
                obj.iNeighbor(iCell_1D) = {vNeighbor};  
                obj.iConnectType(iFaceTemp) = [];     
                iReduce = (obj.iFace >= iFaceTemp);
                obj.iFace(iReduce) = obj.iFace(iReduce) - 1;
            end
        end
        
        %% Check connectivity
        function [isConnected, Aface, Ldist1, Ldist2] = checkConnectivity...
                                  (obj, frac, iCellFrac, grid, ID_MCell_1, ID_MCell_2)
            vert1 = grid.vecVert{ID_MCell_1};
            vert2 = grid.vecVert{ID_MCell_2};
            c1 = [grid.xCell(ID_MCell_1), grid.yCell(ID_MCell_1), grid.zCell(ID_MCell_1)];
            c2 = [grid.xCell(ID_MCell_2), grid.yCell(ID_MCell_2), grid.zCell(ID_MCell_2)];
            vertf = frac.vecCenterLineFrac{iCellFrac};
            vertConnected = zeros(1,3);
            nConVert = 0;
            isConnected = false;
            isFrac = false;
            
            for i=1:length(vert1)
                for j=1:length(vert2)
                    if isequal(vert1(i,:), vert2(j,:))
                        nConVert = nConVert + 1;
                        vertConnected(nConVert, :) = vert1(i, :);
                        if isequal(vert1(i,:), vertf(1,:)) || isequal(vert1(i,:), vertf(2,:))
                            isFrac = true;
                        end
                        break;
                    end
                end
            end
            
            if nConVert == 2     
                isConnected = true;
                v1 = vertConnected(1, :);
                v2 = vertConnected(2, :);
                Aface = norm(v1-v2)*grid.dz(ID_MCell_1);
                [Ldist1, Ldist2] = obj.calcTwoDistance(c1, c2, v1, v2);
            elseif nConVert == 1 && ~isFrac
                isConnected = true;
                v1 = vertConnected(1, :);
                for i = 1:2
                    vf = vertf(i, :);
                    L1f = v1 - vf;
                    for j =1:size(vert2,1)
                        v2 = vert2(j, :);
                        L2f = v2 - vf;
                        if ~isequal(v1, v2) && norm(cross(L1f, L2f)) < 1e-10 % splitted in one matrices                            
                            Aface = norm(L1f)*grid.dz(ID_MCell_1);
                            [Ldist1, Ldist2] = obj.calcTwoDistance(c1, c2, v1, vf);
                            break;
                        end
                    end
                end                
            else
                Aface = 0;
                Ldist1 = 0;
                Ldist2 = 0;                
            end            
        end
        
        %% calculate two distances
        function [L1, L2] = calcTwoDistance(obj, c1, c2, v1, v2)
            vTemp = [1,1];
            v1 = v1(1:2); v2 = v2(1:2);
            c1 = c1(1:2); c2 = c2(1:2);
            vHor = (v1-v2)/norm(v1-v2);
            vVert = vTemp - (vTemp*vHor')*vHor;
            vVert = vVert / norm(vVert);
            
            L1 = 2*abs((c1-v1)*vVert')/norm(vVert);
            L2 = 2*abs((c2-v1)*vVert')/norm(vVert);
        end
        
        %% Delete inactive connections
        function obj = deleteSplitConnection(obj, frac, iCellFrac, grid,...
                                                ID_MCell_1, ID_MCell_2)
            vec_ID_MCell = [ID_MCell_1, ID_MCell_2];
            for i=1:2
                ID_MCell = vec_ID_MCell(i);
                vNeighbor = obj.iNeighbor{ID_MCell};
                for j=1:length(vNeighbor)
                    ID_neighbor = vNeighbor(j);
                    [isConnected, Aface, Ldist1, Ldist2] = obj.checkConnectivity(frac,...
                                                       iCellFrac, grid, ID_MCell, ID_neighbor);

                    if ~isConnected         
                        obj = obj.deleteConnection(ID_MCell, ID_neighbor);
                    else
                        iCell_1D = ID_MCell;
                        iCellNext_1D = ID_neighbor;
                        iFaceTemp = obj.iFace(iCell_1D, iCellNext_1D);
                        obj.areaFace(iFaceTemp) = Aface;
                        obj.distFace(:, iFaceTemp) = [Ldist1, Ldist2];
                        iCell_1D = ID_neighbor;
                        iCellNext_1D = ID_MCell;
                        iFaceTemp = obj.iFace(iCell_1D, iCellNext_1D);
                        obj.areaFace(iFaceTemp) = Aface;
                        obj.distFace(:, iFaceTemp) = [Ldist2, Ldist1];
                    end
                end
            end
        end
        
        %% Make a new connection between splited matrices.
        function obj = splitConnection(obj, grid, frac, iCellFrac, ID_MCell_1, ID_MCell_2) 
            % Shift 
            obj = obj.shiftArray(ID_MCell_1, ID_MCell_2);            
            
            % Copy connection of M1 to M2.
            obj = obj.copyConnection(grid, ID_MCell_1, ID_MCell_2);
                    
            % Delete connections which are separated by M2.
            obj = obj.deleteSplitConnection(frac, iCellFrac, grid,...
                                                ID_MCell_1, ID_MCell_2);
                        
            % Build vacant matrices.
            obj.transGeo = zeros(1, obj.nFace);
            obj.transTot = zeros(2, obj.nFace);
            obj.avgDensity = zeros(2, obj.nFace);
            obj.derivTransTot = zeros(8, obj.nFace);
            obj.derivAvgDensity = zeros(4, obj.nFace);
        end
        
        %% Initialize geometric transmissibility
        function arg = initializeTransGeo(obj)
            % Harmonic average of K
            k = (obj.distFace(1,:)+obj.distFace(2,:))./...
                (obj.distFace(1,:)./obj.permFace(1,:)+obj.distFace(2,:)./obj.permFace(2,:));
            % Distance between 2 cells
            d = (obj.distFace(1,:)+obj.distFace(2,:))/2;
            % Geometric transmissibility
            geo = obj.conv_alpha*obj.areaFace.*(k./d);
            arg = geo;
        end  

        %% Make a connection list for each fracture cell.
        function obj = makeConnection(obj, grid, frac) 
            nFracCell = frac.nTotCell;
            % Build M-F connections
            for iCellFrac = 1:nFracCell
                ifrac = frac.vecIDFrac(1,iCellFrac);
                fType = frac.fracType(ifrac);
                switch fType
                    case 0 % dpdk
                        obj = obj.makeDPDKMFConnection(iCellFrac, grid, frac);
                    case 1 % edfm
                        obj = obj.makeEDFMMFConnection(iCellFrac, grid, frac);
                    otherwise %dfm
                        obj = obj.makeDFMMFConnection(iCellFrac, grid, frac);
                end
            end
            
            % Build F-F connections
            for iCellFrac = 1:nFracCell
                ifrac = frac.vecIDFrac(1,iCellFrac);
                fType = frac.fracType(ifrac);
                switch fType
                    case 0 % dpdk
                        obj = obj.makeDPDKFFConnection(iCellFrac, grid, frac);
                    case 1 % edfm
                        obj = obj.makeEDFMFFConnection(iCellFrac, grid, frac);
                    otherwise %dfm
                        obj = obj.makeDFMFFConnection(iCellFrac, grid, frac);
                end
                obj = obj.makeDifferentFFConnection(iCellFrac, grid, frac);
            end
        end
        
        function obj = makeDFMMFConnection(obj, iCellFrac, grid, frac)
            vID_MCell = frac.vecIDCell{1,iCellFrac};
            ID_FCell = iCellFrac + grid.nMatrix; %obj.vGridIDToFCell(1,ID_MCell);
            ifrac = frac.vecIDFrac(1,iCellFrac);            
            cf = frac.center(:, ifrac)';
            nn = frac.normals(:, ifrac);
            oneAper = frac.aperture(ifrac);
            
            for i=1:length(vID_MCell)
                % Geometric part = k*A/d
                % Distance = norm(center of matrix - center of fracture)
                ID_MCell = vID_MCell(i);
                cm = [grid.xCell(1,ID_MCell), grid.yCell(1,ID_MCell), grid.zCell(1,ID_MCell)];
                assumedDist = 2*abs((cf-cm)*nn); 

                % Area = Lfrac*dz
                Lp_1 = frac.vecCenterLineFrac{iCellFrac}(1,:);
                Lp_2 = frac.vecCenterLineFrac{iCellFrac}(2,:);
                assumedArea = norm(Lp_1-Lp_2)*grid.dz(1,ID_MCell);

                % matrix to frac
                obj.nFace = obj.nFace+ 1;
                obj.iFace(ID_MCell, ID_FCell) = obj.nFace;
                obj.iCell(obj.nFace) = ID_MCell;
                obj.iCellNext(obj.nFace) = ID_FCell;
                obj.iGravity(obj.nFace) = 0;
                obj.areaFace(obj.nFace) = assumedArea;
                obj.distFace(:, obj.nFace) = [assumedDist-oneAper, oneAper];
                % Assume isotropic perm for both matrices and fractures
                obj.permFace(:, obj.nFace) = [grid.kx(ID_MCell), grid.kx(ID_FCell)];
                obj.iNeighbor(ID_MCell) = {[obj.iNeighbor{ID_MCell}, ID_FCell]};
                obj.iConnectType(obj.nFace) = 6;   
                
                % frac to matrix
                obj.nFace = obj.nFace+ 1;
                obj.iFace(ID_FCell, ID_MCell) = obj.nFace;
                obj.iCell(obj.nFace) = ID_FCell;
                obj.iCellNext(obj.nFace) = ID_MCell;
                obj.iGravity(obj.nFace) = 0;
                obj.areaFace(obj.nFace) = assumedArea;
                obj.distFace(:, obj.nFace) = [oneAper, assumedDist-oneAper];
                % Assume isotropic perm for both matrices and fractures
                obj.permFace(:, obj.nFace) = [grid.kx(ID_FCell), grid.kx(ID_MCell)];
                if length(obj.iNeighbor) < ID_FCell
                    obj.iNeighbor(ID_FCell) = {[ID_MCell]};
                else
                    obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell}, ID_MCell]};                    
                end
                obj.iConnectType(obj.nFace) = 6;    
            end
        end
        
        function obj = makeEDFMMFConnection(obj, iCellFrac, grid, frac)
            vecID_MCell = frac.vecIDCell{iCellFrac};
            ID_FCell = iCellFrac + grid.nMatrix; 
            for i=1:length(vecID_MCell)
                ID_MCell = vecID_MCell(i);
                % Geometric part = k*A/d
                assumedDist = grid.calculateDistance(iCellFrac, frac, ID_MCell);

                % Area = Lfrac*dz
                Lp_1 = frac.vecCenterLineFrac{ID_FCell-grid.nMatrix}(1,:);
                Lp_2 = frac.vecCenterLineFrac{ID_FCell-grid.nMatrix}(2,:);
                assumedArea = 2*norm(Lp_1-Lp_2)*grid.dz(1,ID_MCell);
                % matrix to frac
                obj.nFace = obj.nFace+ 1;
                obj.iFace(ID_MCell, ID_FCell) = obj.nFace;
                obj.iCell(obj.nFace) = ID_MCell;
                obj.iCellNext(obj.nFace) = ID_FCell;
                obj.iGravity(obj.nFace) = 0;
                obj.areaFace(obj.nFace) = assumedArea;
                obj.distFace(:, obj.nFace) = [assumedDist, assumedDist];
                % Assume isotropic perm for both matrices and fractures
                obj.permFace(:, obj.nFace) = [grid.kx(ID_MCell), grid.kx(ID_MCell)];
                obj.iNeighbor(ID_MCell) = {[obj.iNeighbor{ID_MCell},ID_FCell]};
                obj.iConnectType(obj.nFace) = 1;            

                % frac to matrix
                obj.nFace = obj.nFace + 1;
                obj.iFace(ID_FCell, ID_MCell) = obj.nFace;
                obj.iCell(obj.nFace) = ID_FCell;
                obj.iCellNext(obj.nFace) = ID_MCell;
                obj.iGravity(obj.nFace) = 0;
                obj.areaFace(obj.nFace) = assumedArea;
                obj.distFace(:, obj.nFace) = [assumedDist, assumedDist];
                % Assume isotropic perm for both matrices and fractures
                obj.permFace(:, obj.nFace) = [grid.kx(ID_MCell), grid.kx(ID_MCell)];
                obj.iNeighbor(ID_FCell) = {[ID_MCell]};
                obj.iConnectType(obj.nFace) = 3;   

                isPEDFM = false;
                if isPEDFM                
                    ifrac = frac.vecIDFrac(iCellFrac);                
                    % find a neighbor
                    LineEDFM = [Lp_1;Lp_2];
                    [ID_MCell_neighbor,assumedDist_n] = obj.findNeighbor(grid, ID_MCell, LineEDFM);
                    obj = obj.deleteConnection(ID_MCell, ID_MCell_neighbor); 

                    % matrix to frac
                    obj.nFace = obj.nFace+ 1;
                    obj.iFace(ID_MCell_neighbor, ID_FCell) = obj.nFace;
                    obj.iCell(obj.nFace) = ID_MCell_neighbor;
                    obj.iCellNext(obj.nFace) = ID_FCell;
                    obj.iGravity(obj.nFace) = 0;
                    obj.areaFace(obj.nFace) = assumedArea/2;
                    obj.distFace(:, obj.nFace) = [2*assumedDist_n-frac.aperture(ifrac), frac.aperture(ifrac)];
                    % Assume isotropic perm for both matrices and fractures
                    obj.permFace(:, obj.nFace) = [grid.kx(ID_MCell_neighbor), grid.kx(ID_FCell)];
                    obj.iNeighbor(ID_MCell_neighbor) = {[obj.iNeighbor{ID_MCell_neighbor},ID_FCell]};
                    obj.iConnectType(obj.nFace) = 1;            

                    % frac to matrix
                    obj.nFace = obj.nFace+ 1;
                    obj.iFace(ID_FCell, ID_MCell_neighbor) = obj.nFace;
                    obj.iCell(obj.nFace) = ID_FCell;
                    obj.iCellNext(obj.nFace) = ID_MCell_neighbor;
                    obj.iGravity(obj.nFace) = 0;
                    obj.areaFace(obj.nFace) = assumedArea/2;
                    obj.distFace(:, obj.nFace) = [frac.aperture(ifrac), 2*assumedDist_n-frac.aperture(ifrac)];
                    % Assume isotropic perm for both matrices and fractures
                    obj.permFace(:, obj.nFace) = [grid.kx(ID_FCell), grid.kx(ID_MCell_neighbor)];
                    obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell},ID_MCell_neighbor]};
                    obj.iConnectType(obj.nFace) = 3;   

                    %grid.vGridIDToMultiFrac(ID_MCell_neighbor) = {[grid.vGridIDToMultiFrac{ID_MCell_neighbor}, ifrac]};
                    %grid.vGridIDToMultiFCell(ID_MCell_neighbor) = {[grid.vGridIDToMultiFCell{ID_MCell_neighbor}, ID_FCell]};
                    %grid.nEmbeddedFrac(ID_MCell_neighbor) = grid.nEmbeddedFrac(ID_MCell_neighbor) + 1; 
                    %grid.vFracType(ID_MCell_neighbor) = {[grid.vFracType{ID_MCell_neighbor}, 2]};
                    %frac.vecIDCell{ID_FCell-grid.nMatrix} = [frac.vecIDCell{ID_FCell-grid.nMatrix}, ID_MCell_neighbor];
                end
            end
        end
        
        
        function obj = makeDFMFFConnection(obj, iCellFrac, grid, frac)
            ID_MCell = frac.vecIDCell{iCellFrac}; 
            ID_FCell = iCellFrac + grid.nMatrix; %obj.vGridIDToFCell(1,ID_MCell);
            vec_M_Neighbor = [obj.iNeighbor{ID_MCell}];
            %vec_M_Neighbor = vec_M_Neighbor(1:(length(vec_M_Neighbor)-grid.nEmbeddedFrac(1,ID_MCell)));
            vec_M_Neighbor(vec_M_Neighbor>grid.nMatrix) = [];
            vec_M_Neighbor = unique(vec_M_Neighbor);
            ID_F_Neighbor = [];
            ID_M_Neighbor = [];
            for i=1:length(vec_M_Neighbor)
                if ~isempty(grid.vGridIDToMultiFCell{vec_M_Neighbor(i)})
                    vF = grid.vGridIDToMultiFCell{vec_M_Neighbor(i)};
                    for j=1:length(vF)
                       if sum(ismember(ID_F_Neighbor, vF(j)))==0
                            ID_M_Neighbor = [ID_M_Neighbor, vec_M_Neighbor(i)];
                            ID_F_Neighbor = [ID_F_Neighbor, vF(j)];
                       end
                    end
                end
            end
            ii = 0;
            while ii < length(ID_F_Neighbor)
                ii = ii + 1;
                if ID_F_Neighbor(ii) <= grid.nMatrix
                    ID_F_Neighbor(ii) = [];
                    ID_M_Neighbor(ii) = [];
                    ii = ii -1;
                end
            end
            ifrac = frac.vecIDFrac(1,iCellFrac);
            
            for i=1:length(ID_F_Neighbor)                
                ifrac_neighbor = frac.vecIDFrac(ID_F_Neighbor(i)-grid.nMatrix);
                if ID_F_Neighbor(i) > 0 && ifrac == ifrac_neighbor% connect ID_FCELL to ID_F_Neighbor(i)
                    
                    % ID_FCELL to ID_F_Neighbor(i)
                    obj.nFace = obj.nFace+ 1;
                    iFaceMM = obj.iFace(ID_MCell, ID_M_Neighbor(i));
                    ID_MCell_temp = ID_MCell(iFaceMM>0);
                    iFaceMM_temp = iFaceMM(iFaceMM>0);
                    ID_MCell_temp = ID_MCell_temp(1); iFaceMM_temp = iFaceMM_temp(1);
                    obj.iFace(ID_FCell, ID_F_Neighbor(i)) = obj.nFace;
                    obj.iCell(obj.nFace) = ID_FCell;
                    obj.iCellNext(obj.nFace) = ID_F_Neighbor(i);
                    obj.iGravity(obj.nFace) = obj.iGravity(iFaceMM_temp);
                    obj.areaFace(obj.nFace) = grid.dz(ID_MCell_temp)*frac.aperture(ifrac);
                    
                    % Distance = norm(center of matrix - center of fracture)
                    c_f1 = [grid.xCell(1,ID_FCell);grid.yCell(1,ID_FCell);grid.zCell(1,ID_FCell)];
                    c_f2 = [grid.xCell(1,ID_F_Neighbor(i));grid.yCell(1,ID_F_Neighbor(i));...
                            grid.zCell(1,ID_F_Neighbor(i))];                        
                    Lp_1 = frac.vecCenterLineFrac{iCellFrac}(1,:);
                    Lp_2 = frac.vecCenterLineFrac{iCellFrac}(2,:);
                    d1 = norm(Lp_1-Lp_2);
                    d2 = 2*norm(c_f1-c_f2) - d1;
                    obj.distFace(:, obj.nFace) = [d1, d2];
                    
                    obj.permFace(:, obj.nFace) = [grid.kx(ID_FCell), grid.kx(ID_F_Neighbor(i))];
                    obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell}, ID_F_Neighbor(i)]};
                    obj.iConnectType(obj.nFace) = 4;   
                end
            end            
        end
        
        
        function obj = makeEDFMFFConnection(obj, iCellFrac, grid, frac)
            ID_MCell = frac.vecIDCell{iCellFrac};   
            ID_MCell = ID_MCell(1);
            ID_FCell = iCellFrac + grid.nMatrix; %obj.vGridIDToFCell(1,ID_MCell);
            vec_M_Neighbor = obj.iNeighbor{ID_MCell};
            vec_M_Neighbor = vec_M_Neighbor(1:(length(vec_M_Neighbor)-grid.nEmbeddedFrac(1,ID_MCell)));
            ID_F_Neighbor = [];
            ID_M_Neighbor = [];
            for i=1:length(vec_M_Neighbor)
                if ~isempty(grid.vGridIDToMultiFCell{vec_M_Neighbor(i)})
                    vF = grid.vGridIDToMultiFCell{vec_M_Neighbor(i)};
                    ID_M_Neighbor = [ID_M_Neighbor, vec_M_Neighbor(i)*ones(1,length(vF))];
                    ID_F_Neighbor = [ID_F_Neighbor, vF];
                end
            end
            
            ii = 0;
            while ii < length(ID_F_Neighbor)
                ii = ii + 1;
                if ID_F_Neighbor(ii) <= grid.nMatrix
                    ID_F_Neighbor(ii) = [];
                    ID_M_Neighbor(ii) = [];
                    ii = ii -1;
                end
            end
            ifrac = frac.vecIDFrac(iCellFrac);
            
            for i=1:length(ID_F_Neighbor)
                ifrac_neighbor = frac.vecIDFrac(ID_F_Neighbor(i)-grid.nMatrix);
                if ID_F_Neighbor(i) > 0 && ifrac == ifrac_neighbor % connect ID_FCELL to ID_F_Neighbor(i)
                    
                    % ID_FCELL to ID_F_Neighbor(i)
                    obj.nFace = obj.nFace+ 1;
                    iFaceMM = obj.iFace(ID_MCell, ID_M_Neighbor(i));
                    obj.iFace(ID_FCell, ID_F_Neighbor(i)) = obj.nFace;
                    obj.iCell(obj.nFace) = ID_FCell;
                    obj.iCellNext(obj.nFace) = ID_F_Neighbor(i);
                    obj.iGravity(obj.nFace) = obj.iGravity(iFaceMM);
                    obj.areaFace(obj.nFace) = grid.dz(ID_MCell)*frac.aperture(ifrac);
                    
                    % Distance = norm(center of matrix - center of fracture)                        
                    c_f1 = [grid.xCell(1,ID_FCell);grid.yCell(1,ID_FCell);grid.zCell(1,ID_FCell)];
                    c_f2 = [grid.xCell(1,ID_F_Neighbor(i));grid.yCell(1,ID_F_Neighbor(i));...
                            grid.zCell(1,ID_F_Neighbor(i))];                        
                    Lp_1 = frac.vecCenterLineFrac{iCellFrac}(1,:);
                    Lp_2 = frac.vecCenterLineFrac{iCellFrac}(2,:);
                    d1 = norm(Lp_1-Lp_2);
                    d2 = 2*norm(c_f1-c_f2)-d1;                    
                    obj.distFace(:, obj.nFace) = [d1, d2];                    
                    
                    obj.permFace(:, obj.nFace) = [grid.kx(ID_FCell), grid.kx(ID_F_Neighbor(i))];
                    obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell}, ID_F_Neighbor(i)]};
                    obj.iConnectType(obj.nFace) = 4;   
                end
            end            
        end
        
        function [ID_MCell_neighbor, assumedDist_n] = findNeighbor(obj, grid, iM, LineEDFM)
            vNeighbor = obj.iNeighbor{iM};     
            L1 = LineEDFM(1,:);
            L2 = LineEDFM(2,:);
            cf = (L1+L2)/2;
            ID_MCell_neighbor = 0;
            assumedDist_n = 0;
            basisLineFrac = L2-L1;
            basisLineFrac = [basisLineFrac(2), -basisLineFrac(1), 0];
            cm = [grid.xCell(1,iM), grid.yCell(1,iM), grid.zCell(1,iM)];
            for i=1:length(vNeighbor)
                iN = vNeighbor(i);
                cn = [grid.xCell(1,iN), grid.yCell(1,iN), grid.zCell(1,iN)];
                basisLineMat = cm-cn;
                basisLineMat = [basisLineMat(2), -basisLineMat(1), 0];
                pointM_to_linefrac = (cm-L1)*basisLineFrac';
                pointN_to_linefrac = (cn-L1)*basisLineFrac';
                pointF1_to_linemat = (L1-cm)*basisLineMat';
                pointF2_to_linemat = (L2-cm)*basisLineMat';
                
                if pointM_to_linefrac*pointN_to_linefrac <= 1e-10 && ...
                        pointF1_to_linemat*pointF2_to_linemat < 1e-10                    
                    assumedDist_n = norm(cf-cn);
                    if abs(pointM_to_linefrac) < 1e-10 && pointN_to_linefrac ~= 0
                        if cm(1) < cn(1) || cm(2) < cn(2)
                            ID_MCell_neighbor = iN;
                            return;
                        end
                    elseif pointM_to_linefrac ~= 0 && pointN_to_linefrac ~= 0
                        ID_MCell_neighbor = iN;
                        return;
                    end
                end
            end
        end
        
        function obj = makeDPDKMFConnection(obj, iCellFrac, grid, frac)
            ID_MCell = frac.vecIDCell{iCellFrac};
            ID_FCell = iCellFrac + grid.nMatrix; %obj.vGridIDToFCell(1,ID_MCell);
            ifrac = frac.vecIDFrac(1,iCellFrac);
            % shape factor
            % kazemi (see KT Lim and Azia, 194), nSetFrac = 2;    % 2d
            shapeFactor = grid.shapeFactor(1,ID_FCell); % dim: 1/L^2
            % Geometric part = k*A/d = k*(ShapeFractor*vol)
            % Assume A = dxdy, d = A/(shape*Vol)=1/(shape*dz)
            assumedArea = grid.dx(1,ID_MCell)*grid.dy(1,ID_MCell)*frac.fracPoro(ifrac);
            assumedDist = 1/(grid.dz(1,ID_MCell)*shapeFactor);
            % matrix to frac
            obj.nFace = obj.nFace+ 1;
            obj.iFace(ID_MCell, ID_FCell) = obj.nFace;
            obj.iCell(obj.nFace) = ID_MCell;
            obj.iCellNext(obj.nFace) = ID_FCell;
            obj.iGravity(obj.nFace) = 0;
            obj.areaFace(obj.nFace) = assumedArea;
            obj.distFace(:, obj.nFace) = [assumedDist, assumedDist];
            % Assume isotropic perm for both matrices and fractures
            obj.permFace(:, obj.nFace) = [grid.kx(ID_MCell), grid.kx(ID_MCell)];
            obj.iNeighbor(ID_MCell) = {[obj.iNeighbor{ID_MCell},ID_FCell]};
            obj.iConnectType(obj.nFace) = 1;            
            
            % frac to matrix
            obj.nFace = obj.nFace+ 1;
            obj.iFace(ID_FCell, ID_MCell) = obj.nFace;
            obj.iCell(obj.nFace) = ID_FCell;
            obj.iCellNext(obj.nFace) = ID_MCell;
            obj.iGravity(obj.nFace) = 0;
            obj.areaFace(obj.nFace) = assumedArea;
            obj.distFace(:, obj.nFace) = [assumedDist, assumedDist];
            % Assume isotropic perm for both matrices and fractures
            obj.permFace(:, obj.nFace) = [grid.kx(ID_FCell), grid.kx(ID_MCell)];
            obj.iNeighbor(ID_FCell) = {[ID_MCell]};
            obj.iConnectType(obj.nFace) = 1;     
        end
        
        function obj = makeDPDKFFConnection(obj, iCellFrac, grid, frac)
            ID_MCell = frac.vecIDCell{iCellFrac};
            ID_FCell = iCellFrac + grid.nMatrix; %obj.vGridIDToFCell(1,ID_MCell);
            vec_M_Neighbor = obj.iNeighbor{ID_MCell};
            ID_M_Neighbor =[];
            ID_F_Neighbor =[];
            for i=1:length(vec_M_Neighbor)
                i_temp = vec_M_Neighbor(i);
                vec_F = grid.vGridIDToMultiFCell{i_temp};
                for j=1:length(vec_F)
                    ID_M_Neighbor =[ID_M_Neighbor, i_temp];
                    ID_F_Neighbor =[ID_F_Neighbor, vec_F(j)];                    
                end
            end
                          
            for i=1:length(ID_F_Neighbor)
                if ID_F_Neighbor(i) > 0 % connect ID_FCELL to ID_F_Neighbor(i)                    
                    % ID_FCELL to ID_F_Neighbor(i)                    
                    ifrac_neighbor = frac.vecIDFrac(1,ID_F_Neighbor(i)-grid.nMatrix);
                    type_neighbor = frac.fracType(ifrac_neighbor);
                    if type_neighbor == 0 % DPDK
                        obj.nFace = obj.nFace+ 1;
                        iFaceMM = obj.iFace(ID_MCell, ID_M_Neighbor(i));
                        obj.iFace(ID_FCell, ID_F_Neighbor(i)) = obj.nFace;
                        obj.iCell(obj.nFace) = ID_FCell;
                        obj.iCellNext(obj.nFace) = ID_F_Neighbor(i);
                        obj.iGravity(obj.nFace) = obj.iGravity(iFaceMM);
                        obj.areaFace(obj.nFace) = obj.areaFace(iFaceMM)*frac.fracPoro(ifrac_neighbor);
                        obj.distFace(:, obj.nFace) = obj.distFace(:, iFaceMM);
                        obj.permFace(:, obj.nFace) = [grid.kx(ID_FCell), grid.kx(ID_F_Neighbor(i))];
                        obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell},ID_F_Neighbor(i)]};
                        obj.iConnectType(obj.nFace) = 2;                             
                    end
                end
            end            
        end
        
        function obj = makeDifferentFFConnection(obj, iCellFrac, grid, frac)
            vecID_MCell = frac.vecIDCell{iCellFrac};
            ID_FCell = iCellFrac + grid.nMatrix; 
            ifrac = frac.vecIDFrac(1,iCellFrac);
            fType = frac.fracType(ifrac);
            ID_FCell_embed = [grid.vGridIDToMultiFCell{vecID_MCell}]; 
            ID_FCell_embed = unique(ID_FCell_embed);
            if ID_FCell_embed(length(ID_FCell_embed)) == ID_FCell
                idtemp = ID_FCell_embed(length(ID_FCell_embed));
                ID_FCell_embed(length(ID_FCell_embed)) = ID_FCell_embed(1);
                ID_FCell_embed(1) = idtemp;
            end
            ifrac_embed = frac.vecIDFrac(1,ID_FCell_embed-grid.nMatrix);  
            fType_embed = frac.fracType(ifrac_embed);             
                
            for i=1:length(ID_FCell_embed)
                if ID_FCell ~= ID_FCell_embed(i)
                    if fType == 0 && fType_embed(i) ~= 0 % DPDK - DFM or EDFM     
                        %for j=1:length(vecID_MCell)                   
                        ID_MCell = frac.vecIDCell{ID_FCell-grid.nMatrix};
                           % ID_MCell = vecID_MCell(j);
                        % Distance = norm(center of matrix - center of fracture)
                        if fType_embed(i) == 1 % EDFM
                            assumedDist1 = grid.calculateDistance(ID_FCell_embed(i)-grid.nMatrix, frac, ID_MCell);
                            assumedDist2 = assumedDist1;
                            k1 = grid.kx(ID_FCell);
                            k2 = grid.kx(ID_FCell);
                        else         
                            cf = frac.center(:, ifrac_embed(i));
                            nn = frac.normals(:, ifrac_embed(i));
                            cm = [grid.xCell(1,ID_MCell); grid.yCell(1,ID_MCell); grid.zCell(1,ID_MCell)];
                            assumedDist1 = 2*abs((cf-cm)'*nn)-frac.aperture(ifrac_embed(i));
                            assumedDist2 = frac.aperture(ifrac_embed(i));
                            k1 = grid.kx(ID_FCell);
                            k2 = grid.kx(ID_FCell_embed(i));
                        end

                        % Area = Lfrac*dz
                        Lp_1 = frac.vecCenterLineFrac{ID_FCell_embed(i)-grid.nMatrix}(1,:);
                        Lp_2 = frac.vecCenterLineFrac{ID_FCell_embed(i)-grid.nMatrix}(2,:);
                        assumedArea = norm(Lp_1-Lp_2)*grid.dz(1,ID_MCell)*frac.fracPoro(ifrac);

                        obj.nFace = obj.nFace+ 1;
                        obj.iFace(ID_FCell, ID_FCell_embed(i)) = obj.nFace;
                        obj.iCell(obj.nFace) = ID_FCell;
                        obj.iCellNext(obj.nFace) = ID_FCell_embed(i);
                        obj.iGravity(obj.nFace) = 0;
                        obj.areaFace(obj.nFace) = assumedArea;
                        obj.distFace(:, obj.nFace) = [assumedDist1, assumedDist2];
                        % Assume isotropic perm for both matrices and fractures
                        obj.permFace(:, obj.nFace) = [k1, k2];
                        obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell}, ID_FCell_embed(i)]};
                        obj.iConnectType(obj.nFace) = 6;   
                        %end
                    elseif fType ~= 0 && fType_embed(i) == 0 % EDFM or DFM - DPDK     
                        %for j=1:length(vecID_MCell)                        
                        ID_MCell = frac.vecIDCell{ID_FCell_embed(i)-grid.nMatrix};%vecID_MCell(j); 
                        if fType == 1 % EDFM    
                            assumedDist1 = grid.calculateDistance(ID_FCell-grid.nMatrix, frac, ID_MCell);
                            assumedDist2 = assumedDist1;
                            k1 = grid.kx(ID_FCell_embed(i));
                            k2 = grid.kx(ID_FCell_embed(i));
                        else % DFM
                            ID_MCell_DPDK = frac.vecIDCell{ID_FCell_embed(i)-grid.nMatrix};
                            cf = frac.center(:, ifrac);
                            nn = frac.normals(:, ifrac);
                            cm = [grid.xCell(1,ID_MCell_DPDK); grid.yCell(1,ID_MCell_DPDK);...
                                  grid.zCell(1,ID_MCell_DPDK)];
                            assumedDist1 = 2*abs((cf-cm)'*nn)-frac.aperture(ifrac);
                            assumedDist2 = frac.aperture(ifrac_embed(i));
                            k1 = grid.kx(ID_FCell_embed(i));
                            k2 = grid.kx(ID_FCell);
                        end

                        % Area = Lfrac*dz
                        Lp_1 = frac.vecCenterLineFrac{ID_FCell-grid.nMatrix}(1,:);
                        Lp_2 = frac.vecCenterLineFrac{ID_FCell-grid.nMatrix}(2,:);
                        assumedArea = norm(Lp_1-Lp_2)*grid.dz(1,ID_MCell)*frac.fracPoro(ifrac_embed(i));

                        obj.nFace = obj.nFace+ 1;
                        obj.iFace(ID_FCell, ID_FCell_embed(i)) = obj.nFace;
                        obj.iCell(obj.nFace) = ID_FCell;
                        obj.iCellNext(obj.nFace) = ID_FCell_embed(i);
                        obj.iGravity(obj.nFace) = 0;
                        obj.areaFace(obj.nFace) = assumedArea;
                        obj.distFace(:, obj.nFace) = [assumedDist1, assumedDist2];
                        % Assume isotropic perm for both matrices and fractures
                        obj.permFace(:, obj.nFace) = [k1, k2];
                        obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell}, ID_FCell_embed(i)]};
                        obj.iConnectType(obj.nFace) = 6;  
                        %end
                    elseif fType ~= 0 && fType_embed(i) ~= 0                         
                        ID_MCell = vecID_MCell(1);                         
                        isDFMEDFM = false;
                        if fType == 2 && fType_embed(i) == 1
                            ID_DFM = ID_FCell;
                            ID_EDFM = ID_FCell_embed(i);
                            ifrac_EDFM = ifrac_embed(i);
                            iCellFrac_DFM = ID_DFM - grid.nMatrix;       
                            
                            ID_MCell_nextEDFM = frac.vecIDCell{iCellFrac_DFM};
                            ID_FCell_nextEDFM = unique([grid.vGridIDToMultiFCell{ID_MCell_nextEDFM}]);
                            ID_FCell_nextEDFM = ID_FCell_nextEDFM(ID_FCell_nextEDFM ~=ID_EDFM);
                            iCellFrac_nextEDFM = ID_FCell_nextEDFM - grid.nMatrix;
                            ifrac_nextEDFM = frac.vecIDFrac(1, iCellFrac_nextEDFM);
                            fType_nextEDFM = frac.fracType(ifrac_nextEDFM);
                            EDFMindex = (ifrac_nextEDFM == ifrac_EDFM); 
                            
                            ID_FCell_nextEDFM = ID_FCell_nextEDFM(EDFMindex);
                            iCellFrac_nextEDFM = iCellFrac_nextEDFM(EDFMindex);
                            ifrac_nextEDFM = ifrac_nextEDFM(EDFMindex);
                            fType_nextEDFM =  fType_nextEDFM(EDFMindex);
                            if length(ID_FCell_nextEDFM) > 0
                                isDFMEDFM = true;
                            end
                        end
                        
                        if isDFMEDFM    % EDFM-DFM   
                            
                            [area_avg, dist_avg, perm_avg,...
                                area_avg2, dist_avg2, perm_avg2] = grid.calcIntersection2...
                                (frac, ID_MCell, ID_FCell, ID_EDFM, ID_FCell_nextEDFM);
                            % DFM
                            obj.nFace = obj.nFace+ 1;
                            obj.iFace(ID_FCell, ID_FCell_embed(i)) = obj.nFace;
                            obj.iCell(obj.nFace) = ID_FCell;
                            obj.iCellNext(obj.nFace) = ID_FCell_embed(i);
                            obj.iGravity(obj.nFace) = 0; 
                            obj.areaFace(obj.nFace) = area_avg;  
                            obj.distFace(:, obj.nFace) = [dist_avg, dist_avg];
                            obj.permFace(:, obj.nFace) = [perm_avg, perm_avg];
                            obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell},ID_FCell_embed(i)]};
                            obj.iConnectType(obj.nFace) = 5; 
                            
                            % EDFM
                            obj.nFace = obj.nFace+ 1;
                            obj.iFace(ID_EDFM, ID_FCell_nextEDFM) = obj.nFace;
                            obj.iCell(obj.nFace) = ID_EDFM;
                            obj.iCellNext(obj.nFace) = ID_FCell_nextEDFM;
                            obj.iGravity(obj.nFace) = 0; 
                            obj.areaFace(obj.nFace) = area_avg2;  
                            obj.distFace(:, obj.nFace) = [dist_avg2, dist_avg2];
                            obj.permFace(:, obj.nFace) = [perm_avg2, perm_avg2];
                            obj.iNeighbor(ID_EDFM) = {[obj.iNeighbor{ID_EDFM},ID_FCell_nextEDFM]};
                            obj.iConnectType(obj.nFace) = 5;
                        elseif fType == 2 && fType_embed(i) == 2 % DFM-DFM,                         
                            obj.nFace = obj.nFace+ 1;
                            obj.iFace(ID_FCell, ID_FCell_embed(i)) = obj.nFace;
                            obj.iCell(obj.nFace) = ID_FCell;
                            obj.iCellNext(obj.nFace) = ID_FCell_embed(i);
                            obj.iGravity(obj.nFace) = 0; 
                            [area_avg, dist_avg, perm_avg] = grid.calculateIntersectionDFM(frac, ID_MCell, ID_FCell, ID_FCell_embed(i));
                            obj.areaFace(obj.nFace) = area_avg;  
                            obj.distFace(:, obj.nFace) = [dist_avg, dist_avg];
                            obj.permFace(:, obj.nFace) = [perm_avg, perm_avg];
                            obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell},ID_FCell_embed(i)]};
                            obj.iConnectType(obj.nFace) = 5;                             
                            
                        else % DFM-half EDFM, EDFM-half DFM                           
                            obj.nFace = obj.nFace+ 1;
                            obj.iFace(ID_FCell, ID_FCell_embed(i)) = obj.nFace;
                            obj.iCell(obj.nFace) = ID_FCell;
                            obj.iCellNext(obj.nFace) = ID_FCell_embed(i);
                            obj.iGravity(obj.nFace) = 0; 
                            [area_avg, dist_avg, perm_avg] = grid.calculateIntersectionDistance(frac, ID_MCell, ID_FCell, ID_FCell_embed(i));
                            obj.areaFace(obj.nFace) = area_avg;  
                            obj.distFace(:, obj.nFace) = [dist_avg, dist_avg];
                            obj.permFace(:, obj.nFace) = [perm_avg, perm_avg];
                            obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell},ID_FCell_embed(i)]};
                            obj.iConnectType(obj.nFace) = 5;  
                        end
                    end
                    
                elseif fType == 2 && fType_embed(i) == 2 && ...
                        length(frac.vecListIntersect{ID_FCell-grid.nMatrix}) > 0% DFM-DFM,                       
                    ID_MCell = vecID_MCell(1);  
                    ID_F_Final = frac.vecListIntersect{ID_FCell-grid.nMatrix} + grid.nMatrix;
                    for j=1:length(ID_FCell_embed)
                        ID_F_Final(ID_F_Final == ID_FCell_embed(j)) = [];
                    end
                    obj.nFace = obj.nFace+ 1;
                    obj.iFace(ID_FCell, ID_F_Final) = obj.nFace;
                    obj.iCell(obj.nFace) = ID_FCell;
                    obj.iCellNext(obj.nFace) = ID_F_Final;
                    obj.iGravity(obj.nFace) = 0; 
                    [area_avg, dist_avg, perm_avg] = grid.calculateIntersectionDFM(frac, ID_MCell, ID_FCell, ID_F_Final);
                    obj.areaFace(obj.nFace) = area_avg;  
                    obj.distFace(:, obj.nFace) = [dist_avg, dist_avg];
                    obj.permFace(:, obj.nFace) = [perm_avg, perm_avg];
                    obj.iNeighbor(ID_FCell) = {[obj.iNeighbor{ID_FCell},ID_F_Final]};
                    obj.iConnectType(obj.nFace) = 5;                                 
                    
                    
                end
            end         
        end
        
        %% Update upwinding index
        function iUpwind = updateUpwind(obj, fluid, grid, avgDensity)            
            % Upwinding considering gravity effect.
            pot1Wat = fluid.Po(obj.iCell) - obj.iGravity.*...
                avgDensity(1,:).*grid.zCell(obj.iCell);
            pot2Wat = fluid.Po(obj.iCellNext) - obj.iGravity.*...
                avgDensity(1,:).*grid.zCell(obj.iCellNext);
            pot1Oil = fluid.Po(obj.iCell) - obj.iGravity.*...
                avgDensity(2,:).*grid.zCell(obj.iCell);
            pot2Oil = fluid.Po(obj.iCellNext) - obj.iGravity.*...
                avgDensity(2,:).*grid.zCell(obj.iCellNext);
            iUpwindWat = pot1Wat>pot2Wat;
            iUpwindOil = pot1Oil>pot2Oil;
            iUpwind = [iUpwindWat;iUpwindOil];
        end                
        %% Update total transmissibility
        function arg = updateTransTot(obj, fluid, iUpwind)
            fluidTermWat = iUpwind(1,:).*(fluid.Krw(obj.iCell)./...
                (fluid.Bw(obj.iCell).*fluid.visW(obj.iCell)))...
                +(1-iUpwind(1,:)).*(fluid.Krw(obj.iCellNext)./...
                (fluid.Bw(obj.iCellNext).*fluid.visW(obj.iCellNext)));
            fluidTermOil = iUpwind(2,:).*(fluid.Kro(obj.iCell)...
                ./(fluid.visO(obj.iCell).*fluid.Bo(obj.iCell)))...
                +(1-iUpwind(2,:)).*(fluid.Kro(obj.iCellNext)...
                ./(fluid.visO(obj.iCellNext).*fluid.Bo(obj.iCellNext)));

            arg = [obj.transGeo.*fluidTermWat; obj.transGeo.*fluidTermOil];            
        end        
        %% Update average density[psi/ft]
        function arg = updateAvgDensity(obj, fluid)
            gamma_denW_avg = (fluid.gamma_denW(obj.iCell)+fluid.gamma_denW(obj.iCellNext))/2;
            gamma_denO_avg = (fluid.gamma_denO(obj.iCell)+fluid.gamma_denO(obj.iCellNext))/2;
            arg = [gamma_denW_avg;gamma_denO_avg];
        end
        
        %% Update derivative of formation volume factor with respect
        %  to Sw and Po of each cell of each face.
        function arg = updateDerivTransTot(obj, fluid, iUpwind)
            %dRw_i/dPo_i
            dWatPo = iUpwind(1,:).*obj.transGeo.*...
                    fluid.Krw(obj.iCell)./fluid.visW(obj.iCell).*(...
                    fluid.deriv_bw(obj.iCell)-fluid.bw(obj.iCell)./...
                        fluid.visW(obj.iCell).*fluid.derivVisW(obj.iCell));                    
            %dRw_i/dPo_next
            dWatPoNext = (1-iUpwind(1,:)).*obj.transGeo.*...
                    fluid.Krw(obj.iCellNext)./fluid.visW(obj.iCellNext).*(...
                    fluid.deriv_bw(obj.iCellNext)-fluid.bw(obj.iCellNext)./...
                        fluid.visW(obj.iCellNext).*fluid.derivVisW(obj.iCellNext));                    
            %dRw_i/dSw_i
            dWatSw = iUpwind(1,:).*obj.transGeo.*fluid.bw(obj.iCell)...
                ./fluid.visW(obj.iCell).*fluid.derivKrw(obj.iCell);            
            %dRw_i/dSw_next
            dWatSwNext = (1-iUpwind(1,:)).*obj.transGeo.*fluid.bw(obj.iCellNext)...
                ./fluid.visW(obj.iCellNext).*fluid.derivKrw(obj.iCellNext);            
            %dRo_i/dPo_i
            dOilPo = iUpwind(2,:).*obj.transGeo.*...
                    fluid.Kro(obj.iCell)./fluid.visO(obj.iCell).*(...
                    fluid.deriv_bo(obj.iCell)-fluid.bo(obj.iCell)./...
                        fluid.visO(obj.iCell).*fluid.derivVisO(obj.iCell));                    
            %dRo_i/dPo_next
            dOilPoNext = (1-iUpwind(2,:)).*obj.transGeo.*...
                    fluid.Kro(obj.iCellNext)./fluid.visO(obj.iCellNext).*(...
                    fluid.deriv_bo(obj.iCellNext)-fluid.bo(obj.iCellNext)./...
                        fluid.visO(obj.iCellNext).*fluid.derivVisO(obj.iCellNext));                     
            %dRw_i/dSw_i
            dOilSw = iUpwind(2,:).*obj.transGeo.*fluid.bo(obj.iCell)...
                ./fluid.visO(obj.iCell).*fluid.derivKro(obj.iCell);            
            %dRw_i/dSw_next
            dOilSwNext = (1-iUpwind(2,:)).*obj.transGeo.*fluid.bo(obj.iCellNext)...
                ./fluid.visO(obj.iCellNext).*fluid.derivKro(obj.iCellNext);
            arg = [dWatPo; dWatPoNext; dWatSw; dWatSwNext;dOilPo; dOilPoNext;...
                dOilSw; dOilSwNext];
        end        
        %% Update derivative of average density
        %[dDenW(iFace)/dPo(iCell), dDenW(iFace)/dPo(iCellNext),
        % dDenO(iFace)/dPo(iCell), dDenO(iFace)/dPo(iCellNext)]
        function arg = updateDerivAvgDensity(obj, fluid)
            [derivDenWat, derivDenOil] = fluid.calcDerivDen(obj.iCell);            
            [derivDenWatNext, derivDenOilNext] = fluid.calcDerivDen(obj.iCellNext);
            arg = [derivDenWat;derivDenWatNext;derivDenOil;derivDenOilNext]/2;
        end        
        %% Update connection properties
        function [iUpwind, transTot, avgDensity, derivTransTot, derivAvgDensity, ...
                  bcTerm, bcDerivTerm]...
                = updateConnnection(obj, fluid, grid)
            avgDensity = updateAvgDensity(obj, fluid);
            derivAvgDensity = updateDerivAvgDensity(obj, fluid);
            iUpwind = updateUpwind(obj, fluid, grid, avgDensity);
            transTot = updateTransTot(obj, fluid, iUpwind);
            derivTransTot = updateDerivTransTot(obj, fluid, iUpwind);
            obj = obj.updateBCterm(grid, fluid);
            bcTerm = obj.bcTerm;
            bcDerivTerm = obj.bcDerivTerm;
        end
        
        function obj = updateBCterm(obj, grid, fluid)
            for iiCell=1:grid.nCell
                if grid.isBoundary(iiCell) == 1 || grid.isBoundary(iiCell)==2
                    Pm = fluid.Po(iiCell);
                    if grid.isBoundary(iiCell) == 1 
                       Pbc = obj.PL;
                       Pup = Pbc;
                    else
                       Pbc = obj.PR;  
                       Pup = Pm;                  
                    end
                    Tgeo = obj.conv_alpha*grid.kx(1,iiCell)*...
                           (grid.dy(1,iiCell)*grid.dz(1,iiCell))/grid.dx(1,iiCell);
                    visO = fluid.visO_std*exp(fluid.compVisOil*(Pup-fluid.Pref));
                    Bo = exp(-fluid.compOil*(Pup-fluid.Pref));
                    Tfluid = 1/(visO*Bo);
                    Ttot = Tgeo*Tfluid;
                    obj.bcTerm(iiCell) = Ttot*(Pbc-Pm);
                    obj.bcDerivTerm(iiCell) = -Ttot;%+Pbc*dTot/dPm;
                    if obj.PL == 0 && obj.PR == 0
                        obj.bcTerm(iiCell) = 0;
                        obj.bcDerivTerm(iiCell) = 0;
                    end
                else
                    obj.bcTerm(iiCell) = 0;
                    obj.bcDerivTerm(iiCell) = 0;
                end
                
                
            end
        end
    end
end

