classdef fracture    
    properties
        center; %ft: 3 by nFrac
        length; %ft: 1 by nFrac
        width; %ft: 1 by nFrac
        strike; %radian: 1 by nFrac, with respect to positive x axis(counterclockwise).
        aperture; %ft: 1 by nFrac
        nFrac;
        fracType; % 1 by nFrac, 0: DPDK, 1:EDFM, 2:DFM
        fracLx; % ft: 1 by nFrac. For DPDK. Describes the x half length of the DPDK domain
        fracLy; % ft: 1 by nFrac. For DPDK. Describes the y half length of the DPDK domain
        fracLz; % ft: 1 by nFrac. For DPDK. Describes the z half length of the DPDK domain
        vecILoc; % i index of the cell with a fracture
        vecJLoc; % j index of the cell with a fracture
        vecKLoc; % k index of the cell with a fracture
        vecCenter; % the center of the cell with a fracture; 3 by # of cells
        vecIDCell; % cell id of the cell with a fracture
        vecIDFrac; % fracture id of the cell with a fracture
        vecNumCell; % 1 by nFrac; number of cells for each fracture
        nTotCell; % the total number of embedded cells
        normals; % 3 by nfrac, normal vectors
        fracPoro; % porosity for DPDK
        fracPerm; % permeability for DPDK
        Lx_matrix; % matrix x length for DPDK
        Ly_matrix; % matrix y length for DPDK
        Lz_matrix; % matrix z length for DPDK
        shapeFactor;
        
        vecCenterLineFrac; % the line passing through fracture center horizontal to the x-y plane
        vecVertFrac;
        vecVolFrac;        
        vecListIntersect;
        vecCenterIntersect;
    end
    
    methods
        function obj = fracture(c, l, w, s, a, n, f, rx, ry, rz, poro, perm, grid)
            obj.center = c;
            obj.length = l;
            obj.width = w;
            obj.strike = s;
            obj.aperture = a;
            obj.nFrac = n;
            obj.fracType = f;            
            obj.fracLx = rx;          
            obj.fracLy = ry;
            obj.fracLz = rz;
            obj.vecNumCell = zeros(1,n);
            obj.nTotCell = 0;
            obj.normals = [sin(obj.strike/180*pi);-cos(obj.strike/180*pi);zeros(1,n)];
            obj.vecCenter = zeros(3,1);
            obj.fracPoro = poro;
            obj.fracPerm = perm;
            obj.vecCenterLineFrac = cell(0);
            obj.vecVertFrac = cell(0);
            obj = obj.definePermPoro(grid);
            obj.vecIDCell = cell(0);
            obj.vecListIntersect = cell(0);
            obj.vecCenterIntersect = cell(0);
        end
        
        function obj = definePermPoro(obj, grid)
            for ifrac = 1:obj.nFrac
                if obj.fracType(ifrac) ~= 0  
                    obj.fracPoro(ifrac) = obj.aperture(ifrac)/obj.length(ifrac); %temporary
                    obj.fracPerm(ifrac) = 1000;%1/12*obj.aperture(ifrac)^2/(1.0623501614639399e-11)*1000; %ft2->md 
                else
                    dx = grid.dx(1);
                    dy = grid.dy(1);
                    dz = grid.dz(1);
                    volFrac = dz*obj.aperture(ifrac)*dx;
                    obj.fracPoro(ifrac) = volFrac/(dx*dy*dz);
                end
            end
        end
        
        function [obj, grid] = initializeFracGeometry(obj, grid)    
            for ifrac = 1:obj.nFrac
                switch obj.fracType(ifrac)
                    case 0 
                        % DPDK: include cells inside the boundary of DPDK domain.
                        % If the cell contains the boundary, exclude the cell.
                        [obj, grid] = obj.initializeFracDPDKGeometry(ifrac, grid);
                    case 1
                        % EDFM: include cells inside the boundary of DPDK domain.
                        % If the cell contains the boundary, exclude the cell.
                        [obj, grid] = obj.initializeFracEDFMGeometry(ifrac, grid);
                    otherwise % DFM                          
                        [obj, grid] = obj.initializeFracDFMGeometry(ifrac, grid);
                end
            end
        end
        % end function
        
        function [obj, grid] = initializeFracEDFMGeometry(obj, ifrac, grid)
            for iCell = 1:grid.nx
                for jCell = 1:grid.ny
                   for kCell = 1:grid.nz
                        iiCell = iCell+(jCell-1)*grid.nx+(kCell-1)*grid.nx*grid.ny;
                        volFrac = 0;
                        [isAbove, centerEDFM, lineEDFM, vertEDFM, volFrac] = obj.calculateCenter(grid, iiCell, ifrac);
                        LfromOrigin = norm(centerEDFM-obj.center(:,ifrac)');
                        if abs(sum(isAbove)) ~= 8 && abs(sum(isAbove)) ~= 0 && LfromOrigin <= obj.length(ifrac)/2
                            obj.nTotCell = obj.nTotCell + 1;
                            obj.vecNumCell(1,ifrac)=obj.vecNumCell(1,ifrac) + 1; 
                            obj.vecILoc(1,obj.nTotCell) = iCell;
                            obj.vecJLoc(1,obj.nTotCell) = jCell;
                            obj.vecKLoc(1,obj.nTotCell) = kCell;
                            obj.vecIDFrac(1,obj.nTotCell) = ifrac;
                            obj.vecCenter(:,obj.nTotCell) = centerEDFM;
                            obj.vecCenterLineFrac(obj.nTotCell) = {lineEDFM};
                            obj.vecVertFrac(obj.nTotCell) = {vertEDFM};
                            obj.vecVolFrac(1,obj.nTotCell) = volFrac;
                            obj.vecIDCell(obj.nTotCell) = {[iiCell]};
                            obj.vecListIntersect(obj.nTotCell) = {[]};
                            obj.vecCenterIntersect(obj.nTotCell) = {[]};
                            if length(grid.vGridIDToMultiFrac) < iiCell                                
                                grid.vGridIDToMultiFrac(iiCell) = {[ifrac]};
                                grid.vGridIDToMultiFCell(iiCell) = {[obj.nTotCell+grid.nMatrix]};                        
                                grid.vFracType(iiCell) = {[obj.fracType(ifrac)]};
                            else
                                grid.vGridIDToMultiFrac(iiCell) = {[grid.vGridIDToMultiFrac{iiCell},ifrac]};
                                grid.vGridIDToMultiFCell(iiCell) = {[grid.vGridIDToMultiFCell{iiCell},...
                                                                     obj.nTotCell+grid.nMatrix]};                        
                                grid.vFracType(iiCell) = {[grid.vFracType{iiCell}, obj.fracType(ifrac)]};
                            end
                        else
                            if length(grid.vGridIDToMultiFrac) < iiCell  
                                grid.vGridIDToMultiFrac(iiCell) = {[]};
                                grid.vGridIDToMultiFCell(iiCell) = {[]};                        
                                grid.vFracType(iiCell) = {[]};
                            end
                        end
                   end
                end
            end                     
        end
        
        function [obj, grid] = splitFracture(obj, grid, ID_MCell_1, ID_MCell_2, iCellFrac)
            iFracGroup = grid.vGridIDToMultiFrac{ID_MCell_1};
            ID_F_Cell_Group = grid.vGridIDToMultiFCell{ID_MCell_1};
            ifrac_DFM = obj.vecIDFrac(1, iCellFrac);
            fType_DFM = obj.fracType(ifrac_DFM);
            ID_FCell_DFM = iCellFrac + grid.nMatrix;
            if fType_DFM ~= 2
                return;
            end
            for i = 1:length(iFracGroup)
                ifrac_split = iFracGroup(i);
                ID_FCell_split = ID_F_Cell_Group(i);
                fType_split = obj.fracType(ifrac_split);
                if fType_split == 0 % DPDK
                    [obj, grid] = obj.addDPDK(grid, ID_MCell_1, ID_MCell_2, ID_FCell_split, ifrac_split); 
                    ID_F_Cell_Group = grid.vGridIDToMultiFCell{ID_MCell_1};            
                elseif fType_split == 1 % EDFM
                    [obj, grid] = obj.addEDFM(grid, ID_MCell_1, ID_MCell_2, ID_FCell_split, ifrac_split);
                    ID_F_Cell_Group = grid.vGridIDToMultiFCell{ID_MCell_1};
                end
            end
        end
        
        function [obj, grid] = addDFM(obj, grid, ID_MCell_1, ID_MCell_2, ID_FCell_split, ifrac_split)   
            vecMCell = [ID_MCell_1, ID_MCell_2];
            vecFCell = [ID_FCell_split, ID_FCell_split+1];
            vecFCellNext = [ID_FCell_split+1, ID_FCell_split];
            [obj, grid] = obj.insertFracture(grid, ID_FCell_split, ifrac_split);
            
            for i=1:2
                iiCell = vecMCell(i);
                ifCell = vecFCell(i);
                ifCellNext = vecFCellNext(i);
                iFracCell = ifCell - grid.nMatrix;
                [isAbove, centerDFM, lineDFM, vertDFM, volFrac] = ...
                    obj.calculateCenter(grid, iiCell, ifrac_split);               
                obj.vecCenter(:,iFracCell) = centerDFM;
                obj.vecCenterLineFrac(iFracCell) = {lineDFM};
                obj.vecVertFrac(iFracCell) = {vertDFM};
                obj.vecVolFrac(1,iFracCell) = volFrac;
                obj.vecIDCell(iFracCell) = {[iiCell]};
                obj.vecListIntersect(iFracCell)  = {[]};   
                obj.vecCenterIntersect(iFracCell)  = {[]};  
                vGrid = grid.vGridIDToMultiFCell{iiCell};
                if sum(ismember(vGrid, ifCell)) == 0
                    vGrid = [vGrid, ifCell];
                    grid.vGridIDToMultiFCell(iiCell) = {vGrid};                    
                end
                
                if sum(ismember(vGrid, ifCellNext)) > 0
                    increaseIndex = (vGrid==ifCellNext);
                    vGrid(increaseIndex) = [];
                    grid.vGridIDToMultiFCell(iiCell) = {vGrid};   
                end
            end  
        end
        
        
        function [obj, grid] = addDPDK(obj, grid, ID_MCell_1, ID_MCell_2, ID_FCell_split, ifrac_split)   
            vecMCell = [ID_MCell_1, ID_MCell_2];
            vecFCell = [ID_FCell_split, ID_FCell_split+1];
            vecFCellNext = [ID_FCell_split+1, ID_FCell_split];
            [obj, grid] = obj.insertFracture(grid, ID_FCell_split, ifrac_split);
            
            for i=1:2
                iiCell = vecMCell(i);
                ifCell = vecFCell(i);
                ifCellNext = vecFCellNext(i);
                iFracCell = ifCell - grid.nMatrix;               
                         
                obj.vecCenter(:,iFracCell) = [grid.xCell(1,iiCell);...
                             grid.yCell(1,iiCell);grid.zCell(1,iiCell)];
                obj.vecCenterLineFrac(iFracCell) = {[]};
                obj.vecVertFrac(iFracCell) = {grid.vecVert{iiCell}};
                obj.vecVolFrac(1,iFracCell) = grid.volCell(1,iiCell);
                obj.vecIDCell(iFracCell) = {[iiCell]};   
                vGrid = grid.vGridIDToMultiFCell{iiCell};
                obj.vecListIntersect(iFracCell)  = {[]};  
                obj.vecCenterIntersect(iFracCell)  = {[]};  
                if sum(ismember(vGrid, ifCell)) == 0
                    vGrid = [vGrid, ifCell];
                    grid.vGridIDToMultiFCell(iiCell) = {vGrid};                    
                end
                
                if sum(ismember(vGrid, ifCellNext)) > 0
                    increaseIndex = (vGrid==ifCellNext);
                    vGrid(increaseIndex) = [];
                    grid.vGridIDToMultiFCell(iiCell) = {vGrid};   
                end
            end  
        end
        
        function [obj, grid] = insertFracture(obj, grid, ID_FCell_split, ifrac_split)   
            obj.vecNumCell(1,ifrac_split)=obj.vecNumCell(1,ifrac_split) + 1;   
            iCellFrac_split = ID_FCell_split-grid.nMatrix;
            iprev = 1:iCellFrac_split;
            inext = (iCellFrac_split+1):obj.nTotCell;  
            obj.nTotCell = obj.nTotCell + 1;
            obj.vecILoc = [obj.vecILoc(iprev), obj.vecILoc(iCellFrac_split),...
                           obj.vecILoc(inext)];
            obj.vecJLoc = [obj.vecJLoc(iprev), obj.vecJLoc(iCellFrac_split),...
                           obj.vecJLoc(inext)];
            obj.vecKLoc = [obj.vecKLoc(iprev), obj.vecKLoc(iCellFrac_split),...
                           obj.vecKLoc(inext)];
            obj.vecIDFrac = [obj.vecIDFrac(iprev), obj.vecIDFrac(iCellFrac_split),...
                             obj.vecIDFrac(inext)];
            obj.vecCenter = [obj.vecCenter(:, iprev), zeros(3,1), obj.vecCenter(:, inext)];
            obj.vecCenterLineFrac = [obj.vecCenterLineFrac(iprev), {[]}, obj.vecCenterLineFrac(inext)];
            obj.vecVertFrac = [obj.vecVertFrac(iprev), {[]}, obj.vecVertFrac(inext)];
            obj.vecVolFrac = [obj.vecVolFrac(iprev), 0, obj.vecVolFrac(inext)];
            obj.vecIDCell = [obj.vecIDCell(iprev), {[]}, obj.vecIDCell(inext)]; 
            obj.vecListIntersect = [obj.vecListIntersect(iprev), {[]}, obj.vecListIntersect(inext)];  
            obj.vecCenterIntersect = [obj.vecCenterIntersect(iprev), {[]}, obj.vecCenterIntersect(inext)];  
            for i=1:length(grid.vGridIDToMultiFCell)
                vGrid = grid.vGridIDToMultiFCell{i};
                increaseIndex = (vGrid>ID_FCell_split);
                vGrid(increaseIndex) = vGrid(increaseIndex)+1;
                grid.vGridIDToMultiFCell(i) = {vGrid};
            end 
            
            for i=1:length(obj.vecListIntersect)
                vIntersect = obj.vecListIntersect{i};
                vIntersect(vIntersect > iCellFrac_split) = vIntersect(vIntersect > iCellFrac_split)+1;
                obj.vecListIntersect(i) = {vIntersect};
            end
        end        
        
        
        
        function [obj, grid] = addEDFM(obj, grid, ID_MCell_1, ID_MCell_2, ID_FCell_split, ifrac_split)   
            vecMCell = [ID_MCell_1, ID_MCell_2];
            vecFCell = [ID_FCell_split, ID_FCell_split+1];
            vecFCellNext = [ID_FCell_split+1, ID_FCell_split];
            [obj, grid] = obj.insertFracture(grid, ID_FCell_split, ifrac_split);
            
            for i=1:2
                iiCell = vecMCell(i);
                ifCell = vecFCell(i);
                ifCellNext = vecFCellNext(i);
                iFracCell = ifCell - grid.nMatrix;
                [isAbove, centerEDFM, lineEDFM, vertEDFM, volFrac] = ...
                    obj.calculateCenter(grid, iiCell, ifrac_split);               
                obj.vecCenter(:,iFracCell) = centerEDFM;
                obj.vecCenterLineFrac(iFracCell) = {lineEDFM};
                obj.vecVertFrac(iFracCell) = {vertEDFM};
                obj.vecVolFrac(1,iFracCell) = volFrac;
                obj.vecIDCell(iFracCell) = {[iiCell]};   
                obj.vecListIntersect(iFracCell)  = {[]};
                obj.vecCenterIntersect(iFracCell)  = {[]};
                vGrid = grid.vGridIDToMultiFCell{iiCell};
                if sum(ismember(vGrid, ifCell)) == 0
                    vGrid = [vGrid, ifCell];
                    grid.vGridIDToMultiFCell(iiCell) = {vGrid};                    
                end
                
                if sum(ismember(vGrid, ifCellNext)) > 0
                    increaseIndex = (vGrid==ifCellNext);
                    vGrid(increaseIndex) = [];
                    grid.vGridIDToMultiFCell(iiCell) = {vGrid};   
                end
            end  
        end
        
        function [obj, grid] = initializeFracDFMGeometry(obj, ifrac, grid)
            for iCell = 1:grid.nx
                for jCell = 1:grid.ny
                   for kCell = 1:grid.nz
                        iiCell = iCell+(jCell-1)*grid.nx+(kCell-1)*grid.nx*grid.ny;
                        volFrac = 0;
                        [isAbove, centerDFM, lineDFM, vertDFM, volFrac] = obj.calculateCenter(grid, iiCell, ifrac);
                        LfromOrigin = norm(centerDFM-obj.center(:,ifrac)');
                        if abs(sum(isAbove)) ~= 8 && abs(sum(isAbove)) ~= 0 && LfromOrigin <= obj.length(ifrac)/2
                            obj.nTotCell = obj.nTotCell + 1;
                            obj.vecNumCell(1,ifrac) = obj.vecNumCell(1,ifrac) + 1;
                            obj.vecILoc(1,obj.nTotCell) = iCell;
                            obj.vecJLoc(1,obj.nTotCell) = jCell;
                            obj.vecKLoc(1,obj.nTotCell) = kCell;
                            obj.vecIDFrac(1,obj.nTotCell) = ifrac;
                            obj.vecCenter(:,obj.nTotCell) = centerDFM;
                            obj.vecCenterLineFrac(obj.nTotCell) = {lineDFM};
                            obj.vecVertFrac(obj.nTotCell) = {vertDFM};
                            obj.vecVolFrac(1,obj.nTotCell) = volFrac;
                            obj.vecIDCell(obj.nTotCell) = {[iiCell]};
                            obj.vecListIntersect(obj.nTotCell)  = {[]};
                            obj.vecCenterIntersect(obj.nTotCell)  = {[]};
                            
                            if length(grid.vGridIDToMultiFrac) < iiCell                                 
                                grid.vGridIDToMultiFrac(iiCell) = {[ifrac]};
                                grid.vGridIDToMultiFCell(iiCell) = {[obj.nTotCell+grid.nMatrix]};                        
                                grid.vFracType(iiCell) = {[obj.fracType(ifrac)]};
                            else                                
                                grid.vGridIDToMultiFrac(iiCell) = {[grid.vGridIDToMultiFrac{iiCell},ifrac]};
                                grid.vGridIDToMultiFCell(iiCell) = {[grid.vGridIDToMultiFCell{iiCell},...
                                                                 obj.nTotCell+grid.nMatrix]};                           
                                grid.vFracType(iiCell) = {[grid.vFracType{iiCell}, obj.fracType(ifrac)]};
                            end                              
                       else
                            if length(grid.vGridIDToMultiFrac) < iiCell  
                                grid.vGridIDToMultiFrac(iiCell) = {[]};
                                grid.vGridIDToMultiFCell(iiCell) = {[]};    
                                grid.vFracType(iiCell) = {[]};
                            end
                       end
                   end
                end
            end                     
        end
        
        function [isAbove, centerEDFM, lineEDFM, vertEDFM, volFrac] = calculateCenter(obj, grid, iiCell, ifrac) 
            volFrac = 0;
            centerEDFM = zeros(1, 3);
            lineEDFM = zeros(2, 3); % line passing through the center horizontally to the x-y plane.
            vertEDFM = zeros(4, 3);
            dz = grid.dz(1, iiCell);
            vertUp = grid.vecVert{iiCell} + [0, 0, dz/2];
            vertDown = grid.vecVert{iiCell} - [0, 0, dz/2];
            cellVertex =[vertUp;vertDown];
            isAbove = ((cellVertex-obj.center(:,ifrac)')*obj.normals(:,ifrac)>0);
            if abs(sum(isAbove)) ~= 8 && abs(sum(isAbove)) ~= 0
                iLine = 0;
                cf = obj.center(:,ifrac);
                nn = obj.normals(:,ifrac);
                vecVert = grid.vecVert{iiCell};
                for iEdge = 1:4 
                    x1 = vecVert(iEdge, :);
                    if iEdge < 4
                        x2 = vecVert(iEdge+1, :);
                    else
                        x2 = vecVert(1, :);
                    end
                    dx12 = x1(1)-x2(1);
                    dy12 = x1(2)-x2(2);
                    innerProd = dx12*nn(1)+dy12*nn(2);
                    if abs(innerProd) > 1e-10
                        con1 = dy12*x2(1)-dx12*x2(2);
                        con2 = nn(1)*cf(1)+nn(2)*cf(2);
                        xCord = (nn(2)*con1+dx12*con2)/innerProd;
                        yCord = (-nn(1)*con1+dy12*con2)/innerProd;
                        zCord = x1(3);
                        if (xCord >= x1(1) && xCord <= x2(1)) || ...
                                (xCord >= x2(1) && xCord <= x1(1))
                            lineEDFM(iLine+1, :) = [xCord, yCord, zCord];
                            if abs(x1(2)-x2(2)) < 1e-10                                
                                halfAper = obj.aperture(ifrac)/2; %/cos(obj.strike(ifrac)/180*pi);               
                                vertEDFM((2*iLine+1):(2*iLine+2), :) = ones(2,1)*lineEDFM(iLine+1, :) + ...
                                                                   [-halfAper,0,0;
                                                                    halfAper,0,0;];               
                            else                           
                                halfAper = obj.aperture(ifrac)/2; %/cos(obj.strike(ifrac)/180*pi);                
                                vertEDFM((2*iLine+1):(2*iLine+2), :) = ones(2,1)*lineEDFM(iLine+1, :) + ...
                                                                   [0,-halfAper,0;
                                                                    0,halfAper,0;]; 
                            end
                            iLine = iLine + 1;
                        end
                    end
                end
                % ld - lu - ru - rd
                vertEDFM = sortrows(vertEDFM);
                tempVert = vertEDFM(4,:);
                vertEDFM(4,:) = vertEDFM(3,:);
                vertEDFM(3,:) = tempVert;
                
                volFrac = dz*obj.aperture(ifrac)*norm(lineEDFM(1, :)-lineEDFM(2, :));
                centerEDFM = (lineEDFM(1,:) + lineEDFM(2,:)) / 2;
            end
        end
        
        function [obj, grid] = initializeFracDPDKGeometry(obj, ifrac, grid)
            for iCell = 1:grid.nx
                for jCell = 1:grid.ny
                    for kCell = 1:grid.nz
                        iiCell = iCell+(jCell-1)*grid.nx+(kCell-1)*grid.nx*grid.ny;
                        xCell = [grid.xCell(1,iiCell);grid.yCell(1,iiCell);grid.zCell(1,iiCell)];
                        distCellFrac = abs(obj.center(:,ifrac)-xCell); 
                        if distCellFrac(1) + grid.dx(1,iiCell)/2 <= obj.fracLx(ifrac) ...
                            && distCellFrac(2) + grid.dy(1,iiCell)/2 <= obj.fracLy(ifrac) ...
                            && distCellFrac(3) + grid.dz(1,iiCell)/2 <= obj.fracLz(ifrac)   
                            dx = grid.dx(1, iiCell);
                            dz = grid.dz(1, iiCell);                         
                            obj.nTotCell = obj.nTotCell + 1;
                            obj.vecNumCell(1,ifrac)=obj.vecNumCell(1,ifrac) + 1; 
                            obj.vecILoc(1,obj.nTotCell) = iCell;
                            obj.vecJLoc(1,obj.nTotCell) = jCell;
                            obj.vecKLoc(1,obj.nTotCell) = kCell;
                            obj.vecIDFrac(1,obj.nTotCell) = ifrac;
                            obj.vecCenter(:,obj.nTotCell) = xCell;
                            obj.vecIDCell(obj.nTotCell) = {[iiCell]};
                            obj.vecCenterLineFrac(obj.nTotCell) = {[]};
                            obj.vecVertFrac(obj.nTotCell) = {[]};
                            volFrac = dz*obj.aperture(ifrac)*dx;
                            obj.vecVolFrac(1,obj.nTotCell) = volFrac;
                            obj.vecListIntersect(obj.nTotCell)  = {[]};
                            obj.vecCenterIntersect(obj.nTotCell)  = {[]};
                            
                            if length(grid.vGridIDToMultiFrac) < iiCell                             
                                grid.vGridIDToMultiFrac(iiCell) = {[ifrac]};
                                grid.vGridIDToMultiFCell(iiCell) = {[obj.nTotCell+grid.nMatrix]};                        
                                grid.vFracType(iiCell) = {[obj.fracType(ifrac)]};
                            else        
                                grid.vGridIDToMultiFrac(iiCell) = {[grid.vGridIDToMultiFrac{iiCell},ifrac]};
                                grid.vGridIDToMultiFCell(iiCell) = {[grid.vGridIDToMultiFCell{iiCell},...
                                obj.nTotCell+grid.nMatrix]};                           
                                grid.vFracType(iiCell) = {[grid.vFracType{iiCell}, obj.fracType(ifrac)]};
                            end  
                        else
                            if length(grid.vGridIDToMultiFrac) < iiCell  
                                grid.vGridIDToMultiFrac(iiCell) = {[]};
                                grid.vGridIDToMultiFCell(iiCell) = {[]};    
                                grid.vFracType(iiCell) = {[]};
                            end
                        end
                    end
                end                        
            end  
        end
                 
        function [obj, grid] = removeFracCell(obj, grid, iCellFracRemove)
            iCellFracRemove = iCellFracRemove - grid.nMatrix;    
            ifrac = obj.vecIDFrac(1, iCellFracRemove);                
            obj.vecILoc(iCellFracRemove) = [];
            obj.vecJLoc(iCellFracRemove) = [];
            obj.vecKLoc(iCellFracRemove) = [];
            obj.vecIDFrac(iCellFracRemove) = [];
            obj.vecCenter(:,iCellFracRemove) = [];
            obj.vecIDCell(iCellFracRemove) = []; 
            obj.vecListIntersect(iCellFracRemove)  = [];
            obj.vecCenterIntersect(iCellFracRemove)  = [];
            if length(obj.vecCenterLineFrac) >= iCellFracRemove
                obj.vecCenterLineFrac(iCellFracRemove) = [];
                obj.vecVertFrac(iCellFracRemove) = [];
                obj.vecVolFrac(iCellFracRemove) = [];
            end
            obj.vecNumCell(1, ifrac)=obj.vecNumCell(1, ifrac) - 1; 
            obj.nTotCell = obj.nTotCell - 1;   
            iCellFracRemove = iCellFracRemove + grid.nMatrix;  
            
            for i=1:length(grid.vGridIDToMultiFrac)
                vGridToFCell = [grid.vGridIDToMultiFCell{i}];
                reduceIndex = (vGridToFCell>iCellFracRemove);
                removeIndex = (vGridToFCell==iCellFracRemove);
                vGridToFCell(reduceIndex) = vGridToFCell(reduceIndex) - 1;
                vGridToFCell(removeIndex) = [];
                grid.vGridIDToMultiFCell(i) = {vGridToFCell};
                
                vGridToFrac = grid.vGridIDToMultiFrac{i};
                vGridToFrac(removeIndex) = [];
                grid.vGridIDToMultiFrac(i) = {vGridToFrac};    
                
                vFracType = grid.vFracType{i};
                vFracType(removeIndex) = [];
                grid.vFracType(i) = {vFracType};
            end 
        end        
    end
end

