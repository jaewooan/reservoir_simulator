classdef grids
    %grid: Stores rock and geometric properties on each cell.
    
    properties
        %% Geometric property
        Lx = 3300; % Reservoir size in x direction [ft]
        Ly = 110; % Reservoir size in y direction [ft]
        Lz = 1650; % Reservoir size in z direction [ft]
        nx = 0; % the number of grids in x direction
        ny = 0; % the number of grids in y direction
        nz = 0; % the number of grids in z direction
        volCell=0; % Bulk volume of each cell[ft^3]
        z_top=6000; % Top location of the top reservoir[ft]
        dx = 0; % Grid size in x direction [ft]
        dy = 0; % Grid size in y direction [ft]
        dz = 0; % Grid size in z direction [ft]
        xCell; % X location of center of a cell[ft]
        yCell; % Y location of center of a cell[ft]
        zCell; % Z location of center of a cell[ft]
        vecVert;
        nMatrix;
        nFrac;
        nCell;
        
        %% Rock property
        kx; % Permeability in x direction [md]
        ky; % Permeability in y direction [md]
        kz; % Permeability in z direction [md]
        poro; % Porosity of each cell at time level n+1
        poroPrev; % Porosity of each cell at time level n
        compRock; % Rock compressibility[psi-1]
        poro_ref; % Reference porosity for updating porosity
        poro_P_ref = 14.7; % Reference pressure for updating porosity[psi]
        derivPoro; % d(Porosity)/d(Po) at (i,j,k) at time n+1
        
        %% Fracture
        vGridIDToMultiFrac;
        vGridIDToMultiFCell;
        nEmbeddedFrac; % 1 by (nx*ny*nz) ;the number of fractures embedded in a cell
        vFracType;
        shapeFactor;
        
        %% pLOT
        mapPlotToGrid;
        mapPermPlotToGrid;
        mapPlotDist;
        mapOldToNew;
        isBoundary;
    end
    
    methods
        function obj = grids(nGridx,nGridy,nGridz)
            obj.volCell = zeros(1, nGridx*nGridy*nGridz);
            obj.kx = 300*ones(1, nGridx*nGridy*nGridz);
            obj.ky = 300*ones(1, nGridx*nGridy*nGridz);
            obj.kz = 300*ones(1, nGridx*nGridy*nGridz);
            obj.poro = zeros(1, nGridx*nGridy*nGridz);
            obj.poroPrev = zeros(1, nGridx*nGridy*nGridz);
            obj.derivPoro = zeros(1, nGridx*nGridy*nGridz);
            obj.nx = nGridx; 
            obj.ny = nGridy;
            obj.nz = nGridz;
            obj.nMatrix = obj.nx*obj.ny*obj.nz;
            obj.nFrac = 0;
            obj.nCell = obj.nMatrix;
            obj.vecVert = cell(0);
            obj.vGridIDToMultiFrac = cell(0);
            obj.vGridIDToMultiFCell = cell(0);
            obj.nEmbeddedFrac = zeros(1, nGridx*nGridy*nGridz);
            obj.vFracType = cell(0);
            obj.shapeFactor = zeros(1, nGridx*nGridy*nGridz);
            obj.mapPlotToGrid = zeros(nGridx,nGridy,nGridz);
            obj.mapPermPlotToGrid = zeros(nGridx,nGridy,nGridz);
            obj.mapPlotDist = zeros(nGridx,nGridy,nGridz);
        end        
        
        %% Initialize grid
        function obj= initializeGrid(obj)           
            obj.Lx = sum(obj.dx(1)*obj.nx);
            obj.Ly = sum(obj.dy(1)*obj.ny);
            obj.Lz = sum(obj.dz(1)*obj.nz);
            obj.volCell= obj.dx.*obj.dy.*obj.dz;
            % Store the location of each cell
            for i=1:obj.nx
                for j=1:obj.ny
                    for k=1:obj.nz      
                        iCell = i+(j-1)*obj.nx+(k-1)*obj.nx*obj.ny;
                        obj.xCell(1,iCell) = obj.dx(1,iCell)*(i-1/2);
                        obj.yCell(1,iCell) = obj.dy(1,iCell)*(j-1/2);
                        obj.zCell(1,iCell) = obj.z_top+obj.dz(1,iCell)*(k-1/2);
                        obj.vecVert(iCell) =... % ld - lu - ru - rd
                            {[obj.dx(1,iCell)*(i-1), obj.dy(1,iCell)*(j-1), obj.zCell(1,iCell);...
                              obj.dx(1,iCell)*(i-1), obj.dy(1,iCell)*j, obj.zCell(1,iCell);...
                              obj.dx(1,iCell)* i, obj.dy(1,iCell)*j, obj.zCell(1,iCell);...
                              obj.dx(1,iCell)* i, obj.dy(1,iCell)*(j-1), obj.zCell(1,iCell)]};
                        obj.mapOldToNew(iCell) = iCell;
                        if i==1
                            obj.isBoundary(1,iCell) = 1;
                        elseif i==obj.nx
                            obj.isBoundary(1,iCell) = 2;                            
                        end
                    end
                end
            end 
        end
        
        %% Mapping index
        function obj = mappingGridtoPlot(obj)
            % define minimum dx, dy, dz
            % calculate total number of grids.
            nx = obj.nx;
            ny = obj.ny;
            nz = obj.nz;            
            dx = obj.dx(1);
            dy = obj.dy(1);
            dz = obj.dz(1);            
            obj.mapPlotToGrid = zeros(nx, ny, nz);
            
            % 
            for i=1:obj.nMatrix
                xc = [obj.xCell(i), obj.yCell(i), obj.zCell(i)];
                ii = floor(obj.xCell(i)/dx)+1;
                jj = floor(obj.yCell(i)/dy)+1;
                kk = floor(obj.zCell(i)/dz)+1;
                xref = [dx*(ii-1/2), dy*(jj-1/2), dz*(kk-1/2)];
                distance = norm(xref-xc);
                
%                
                if distance < obj.mapPlotDist(ii,jj,kk) || obj.mapPlotToGrid(ii,jj,kk)==0
                      obj.mapPlotDist(ii,jj,kk) = distance;
                      obj.mapPlotToGrid(ii,jj,kk) = i;
                      obj.mapPermPlotToGrid(ii,jj,kk) = i;
                end
                if obj.nEmbeddedFrac(1,i) > 0 && obj.mapPermPlotToGrid(ii,jj,kk) <= obj.nMatrix
                    vFrac = obj.vGridIDToMultiFCell{i};
                    obj.mapPermPlotToGrid(ii,jj,kk) = vFrac(1);
                end
                
            end
            
        end
        
        %% Update fracture porosity, permeability, and grid information.
        function [obj, frac, connection] = defineFractureCells(obj, frac, connection)  
            nFracCell = frac.nTotCell;
            iCellFrac = 0;
            while iCellFrac < nFracCell
                % Define a fracture cell in the grid.
                iCellFrac = iCellFrac + 1;
                if iCellFrac == 228
                    aa=1;
                end
                if  length(frac.vecListIntersect{iCellFrac}) == 0   
                    [obj, frac, connection] = obj.defineFracCell(connection, iCellFrac, frac);
                end
                if nFracCell ~= frac.nTotCell
                    try
                    ifrac = frac.vecIDFrac(1, iCellFrac); 
                    catch
                        aa=1;
                    end
                        
                    if frac.fracType(ifrac) == 0
                        iCellFrac = iCellFrac  - 1;
                    end
                    nFracCell = frac.nTotCell;
                end 
            end
        end 
        
        function [obj, frac, connection] = defineFracCell(obj, connection, iCellFrac, frac) 
            ID_MCell = frac.vecIDCell{iCellFrac};   
            ID_FCell = iCellFrac + obj.nMatrix;
            ifrac = frac.vecIDFrac(1, iCellFrac); 
            if frac.fracType(ifrac) == 0  
                repIndex = ismember(obj.vFracType{ID_MCell}, 0);
                vecIfrac_rep = obj.vGridIDToMultiFrac{ID_MCell};
                vecIfrac_rep = vecIfrac_rep(repIndex);
                ifrac_rep = vecIfrac_rep(1); 
                vecFcell = [obj.vGridIDToMultiFCell{ID_MCell}];
                vecFcell = vecFcell(repIndex);
                ID_FCell_rep = vecFcell(1); 
                if ifrac_rep ~= ifrac
                    [obj, frac] = obj.addDPDKtoDPDK(frac, ID_MCell, ID_FCell, ifrac,...
                                                    ID_FCell_rep);
                else 
                    obj = obj.addDPDK(frac, ID_MCell, ID_FCell, ifrac);                    
                end
            elseif frac.fracType(ifrac) == 1  % EDFM
                obj = obj.addEDFM(frac, ID_MCell, ID_FCell, ifrac);
            else  % DFM 
                repIndex = ismember(obj.vFracType{ID_MCell(1)}, 2);           
                vecIfrac_rep = obj.vGridIDToMultiFrac{ID_MCell(1)};
                vecIfrac_rep = vecIfrac_rep(repIndex);
                ifrac_rep = vecIfrac_rep(1); 
                vecFcell = [obj.vGridIDToMultiFCell{ID_MCell(1)}];
                vecFcell = vecFcell(repIndex);
                %minifrac = min(frac.vecIDFrac(vecFcell-obj.nMatrix));
                %vecFcell(frac.fracType(frac.vecIDFrac(vecFcell-obj.nMatrix))~=2)=[];
                %vecFcell(frac.vecIDFrac(vecFcell-obj.nMatrix)~=minifrac)=[];
                ID_FCell_rep = vecFcell(1); 
                %ifrac_rep = minifrac;
                if ifrac ~= ifrac_rep % intersection       
                    [obj, frac, connection] = obj.addDFMtoDFM(frac, connection, ID_FCell_rep, ID_FCell, ifrac);
                else
                    [obj, frac, connection] = obj.addDFM(frac, connection, ID_MCell, ID_FCell, ifrac);
                end
            end
        end         
        
        
        %% Define an EDFM cell
        function obj = addEDFM(obj, frac, ID_MCell, ID_FCell, ifrac)
            iCellFrac = ID_FCell - obj.nMatrix;
            obj.nCell = obj.nCell + 1; % add a fracture cell
            obj.nFrac = obj.nFrac + 1; % add a fracture cell 
            obj.nEmbeddedFrac(1,ID_MCell) = obj.nEmbeddedFrac(1,ID_MCell) + 1;
            obj.poro(1, ID_FCell) = 1;
            obj.volCell(1, ID_FCell) = frac.vecVolFrac(1,ID_FCell-obj.nMatrix);
            obj.kx(1, ID_FCell) = frac.fracPerm(1,ifrac);
            obj.ky(1, ID_FCell) = frac.fracPerm(1,ifrac);
            obj.kz(1, ID_FCell) = frac.fracPerm(1,ifrac);
            obj.poroPrev(1, ID_FCell) = 0;
            obj.derivPoro(1, ID_FCell) = 0;            
            obj.xCell(1,ID_FCell) = frac.vecCenter(1,ID_FCell-obj.nMatrix);
            obj.yCell(1,ID_FCell) = frac.vecCenter(2,ID_FCell-obj.nMatrix);
            obj.zCell(1,ID_FCell) = frac.vecCenter(3,ID_FCell-obj.nMatrix);
            obj.vecVert(ID_FCell) = {frac.vecVertFrac{iCellFrac}};
            obj.nEmbeddedFrac(ID_FCell) = 0;
            obj.vFracType(ID_FCell) = {[]};
            obj.shapeFactor = 0;
            obj.vGridIDToMultiFrac(ID_FCell) = {[]};
            obj.vGridIDToMultiFCell(ID_FCell) = {[]};
            obj.isBoundary(ID_FCell) = 0;
        end
        
        %% Define an DFM cell
        function [obj, frac, connection] = addDFM(obj, frac, connection, ID_MCell, ID_FCell, ifrac)
            iCellFrac = ID_FCell - obj.nMatrix;
            [obj, frac, connection] = obj.splitMatrix(frac, connection, ID_MCell, iCellFrac, ifrac);
            obj.nCell = obj.nCell + 1; % add a fracture cell
            obj.nFrac = obj.nFrac + 1; % add a fracture cell 
            ID_FCell = iCellFrac + obj.nMatrix;
            
            % fracture
            obj.poro(1, ID_FCell) = 1;
            obj.volCell(1, ID_FCell) = frac.vecVolFrac(iCellFrac);            
            obj.dx(1, ID_FCell) = 0;
            obj.dy(1, ID_FCell) = 0;
            obj.dz(1, ID_FCell) = 0;
            obj.kx(1, ID_FCell) = frac.fracPerm(1,ifrac);
            obj.ky(1, ID_FCell) = frac.fracPerm(1,ifrac);
            obj.kz(1, ID_FCell) = frac.fracPerm(1,ifrac);
            obj.poroPrev(1, ID_FCell) = 0;
            obj.derivPoro(1, ID_FCell) = 0;            
            obj.xCell(1,ID_FCell) = frac.vecCenter(1, iCellFrac);
            obj.yCell(1,ID_FCell) = frac.vecCenter(2, iCellFrac);
            obj.zCell(1,ID_FCell) = frac.vecCenter(3, iCellFrac);
            obj.vecVert(ID_FCell) = {frac.vecVertFrac{iCellFrac}};
            obj.nEmbeddedFrac(ID_FCell) =  0;
            obj.vFracType(ID_FCell) = {[]};
            obj.shapeFactor(ID_FCell) = 0;
            obj.vGridIDToMultiFrac(ID_FCell) = {[]};
            obj.vGridIDToMultiFCell(ID_FCell) = {[]};
            obj.isBoundary(ID_FCell) = 0;
        end

        %% Merge additional DFM cells to the first DFM cell.
        function [obj, frac, connection] = addDFMtoDFM(obj, frac, connection, ID_FCell_1, ID_FCell_2, ifrac_2) 
            iCellFrac1 = ID_FCell_1 - obj.nMatrix;
            iCellFrac2 = ID_FCell_2 - obj.nMatrix;
            vTemp = frac.vecIDCell{iCellFrac1};
            ID_MCell_1 = vTemp(1);
            ID_MCell_2 = vTemp(2);
            [obj, frac, connection] = obj.splitMatrix(frac, connection, ID_MCell_1, iCellFrac2, ifrac_2);                         
            vecID = [frac.vecIDCell{iCellFrac1}, ID_MCell_1, ID_MCell_1+1];
            frac.vecIDCell(iCellFrac1) = {unique(vecID)};              
            ID_MCell_2 = ID_MCell_2 + 1;
            ID_FCell_1 = ID_FCell_1 + 1;
            ID_FCell_2 = ID_FCell_2 + 1;
            [obj, frac, connection] = obj.splitMatrix(frac, connection, ID_MCell_2, iCellFrac2, ifrac_2);                   
            vecID = [frac.vecIDCell{iCellFrac1}, ID_MCell_2, ID_MCell_2+1];
            frac.vecIDCell(iCellFrac1) = {unique(vecID)};  
            ID_FCell_1 = ID_FCell_1 + 1;
            ID_FCell_2 = ID_FCell_2 + 1;
            % fracture
            obj.poro(1, ID_FCell_2) = 1;
            obj.volCell(1, ID_FCell_2) = frac.vecVolFrac(iCellFrac2);        
            obj.dx(1, ID_FCell_2) = 0;
            obj.dy(1, ID_FCell_2) = 0;
            obj.dz(1, ID_FCell_2) = 0;
            obj.kx(1, ID_FCell_2) = frac.fracPerm(1,ifrac_2);
            obj.ky(1, ID_FCell_2) = frac.fracPerm(1,ifrac_2);
            obj.kz(1, ID_FCell_2) = frac.fracPerm(1,ifrac_2);
            obj.poroPrev(1, ID_FCell_2) = 0;
            obj.derivPoro(1, ID_FCell_2) = 0;            
            obj.xCell(1,ID_FCell_2) = frac.vecCenter(1, iCellFrac2);
            obj.yCell(1,ID_FCell_2) = frac.vecCenter(2, iCellFrac2);
            obj.zCell(1,ID_FCell_2) = frac.vecCenter(3, iCellFrac2);            
            obj.nEmbeddedFrac(ID_FCell_2) = 0;            
            obj.vecVert(ID_FCell_2) = {frac.vecVertFrac{iCellFrac2}};
            obj.vFracType(ID_FCell_2) = {[]};
            obj.shapeFactor(1, ID_FCell_2) = 0;
            obj.vGridIDToMultiFrac(ID_FCell_2) = {[]};
            obj.vGridIDToMultiFCell(ID_FCell_2) = {[]};
            obj.isBoundary(ID_FCell_2) = 0;
            
            obj.nCell = obj.nCell + 1; % add a fracture cell
            obj.nFrac = obj.nFrac + 1; % add a fracture cell 
            
            [obj, frac, connection] = obj.splitDFM(frac, connection, iCellFrac1, iCellFrac2); 
        end
        
        %% Split a intersection dfm
        function [obj, frac, connection] = splitDFM(obj, frac, connection, iCellFrac1, iCellFrac2)   
            [center, centerFrac, lineDFM, vertDFM, volFrac] = obj.calculateCenter2(frac, iCellFrac1, iCellFrac2);       
            vecIcellFrac = [iCellFrac1, iCellFrac2+1];
            vecIcellFrac_rev = [iCellFrac2, iCellFrac1];
            vecMCell = unique([frac.vecIDCell{[iCellFrac1,iCellFrac2]}]);
            for k = 1:length(vecMCell)               
                vGrid = obj.vGridIDToMultiFCell{vecMCell(k)};
                vGrid(vGrid == iCellFrac1 + obj.nMatrix) = [];
                vGrid(vGrid == iCellFrac2 + obj.nMatrix) = [];
                obj.vGridIDToMultiFCell(vecMCell(k)) = {vGrid};  
            end
            for k = 1:2
                vecMCell = frac.vecIDCell{vecIcellFrac(k)};%[ID_MCell_1, ID_MCell_2];    
                ID_FCell_split = vecIcellFrac(k) + obj.nMatrix;        
                vecFCell = [ID_FCell_split, ID_FCell_split+1];
                vecFCellNext = [ID_FCell_split+1, ID_FCell_split];
                ifrac_split = frac.vecIDFrac(vecIcellFrac(k)); 
                [frac, obj] = frac.insertFracture(obj, ID_FCell_split, ifrac_split); 
                for i=1:2                
                    vecMCell2 = obj.selectContactCell(vecMCell, lineDFM{2*(k-1)+i});
                    ifCell = vecFCell(i);
                    iFracCell = ifCell - obj.nMatrix;        
                    frac.vecCenterIntersect(iFracCell) = {center};  
                    frac.vecCenter(:,iFracCell) = centerFrac(2*(k-1)+i, :)';
                    frac.vecCenterLineFrac(iFracCell) = {lineDFM{2*(k-1)+i}};
                    frac.vecVertFrac(iFracCell) = {vertDFM{2*(k-1)+i}};
                    frac.vecVolFrac(1,iFracCell) = volFrac(2*(k-1)+i);
                    for j=1:length(vecMCell)
                        vecIDCell_temp = frac.vecIDCell{iFracCell};
                        vecIDCell_temp(vecIDCell_temp  == vecMCell(j)) = []; 
                        frac.vecIDCell(iFracCell) = {vecIDCell_temp};
                    end
                    for j=1:length(vecMCell2)
                        iiCell = vecMCell2(j);
                        frac.vecIDCell(iFracCell) = {[frac.vecIDCell{iFracCell}, iiCell]};   
                        vGrid = obj.vGridIDToMultiFCell{iiCell};
                        if sum(ismember(vGrid, ifCell)) == 0
                            vGrid = [vGrid, ifCell];
                            obj.vGridIDToMultiFCell(iiCell) = {vGrid};                    
                        end
                    end
                end
                
                       
                [obj, frac] = obj.insertFracture(frac, ID_FCell_split, ID_FCell_split+1, 0);
                ID_FCell_1 = ID_FCell_split;
                ID_FCell_2 = ID_FCell_split+1;
                
                % Update the matrix 1
                iFracCell1 = ID_FCell_1 - obj.nMatrix;
                obj.volCell(1,ID_FCell_1) = frac.vecVolFrac(1,iFracCell1);
                obj.xCell(1,ID_FCell_1) = frac.vecCenter(1,iFracCell1);
                obj.yCell(1,ID_FCell_1) = frac.vecCenter(2,iFracCell1);
                obj.zCell(1,ID_FCell_1) = frac.vecCenter(3,iFracCell1);
                obj.vecVert(ID_FCell_1) = {frac.vecVertFrac{iFracCell1}};            

                % Update the matrix 2
                iFracCell2 = ID_FCell_2 - obj.nMatrix;
                obj.volCell(1,ID_FCell_2) = frac.vecVolFrac(1,iFracCell2);
                obj.xCell(1,ID_FCell_2) = frac.vecCenter(1,iFracCell2);
                obj.yCell(1,ID_FCell_2) = frac.vecCenter(2,iFracCell2);
                obj.zCell(1,ID_FCell_2) = frac.vecCenter(3,iFracCell2);
                obj.vecVert(ID_FCell_2) = {frac.vecVertFrac{iFracCell2}};
                obj.kx(1,ID_FCell_2) = obj.kx(1,ID_FCell_1);
                obj.ky(1,ID_FCell_2) = obj.ky(1,ID_FCell_1);
                obj.kz(1,ID_FCell_2) = obj.kz(1,ID_FCell_1);
                obj.poro(1,ID_FCell_2) = obj.poro(1,ID_FCell_1);

                obj.vGridIDToMultiFrac(ID_FCell_2) = {[obj.vGridIDToMultiFrac{ID_FCell_1}]};
                obj.vGridIDToMultiFCell(ID_FCell_2) = {[obj.vGridIDToMultiFCell{ID_FCell_1}]};
                obj.nEmbeddedFrac(1,ID_FCell_2) = obj.nEmbeddedFrac(1,ID_FCell_1);
                obj.vFracType(ID_FCell_2) = {[obj.vFracType{ID_FCell_1}]};
                obj.shapeFactor(1,ID_FCell_2) = obj.shapeFactor(1,ID_FCell_1);
            end              
            
            vListIntersect = [iCellFrac1,iCellFrac1+1,iCellFrac2+1,iCellFrac2+2];
            for i=1:length(vListIntersect)
                frac.vecListIntersect(vListIntersect(i))  = {vListIntersect}; 
            end
        end
        
        function iiCell = selectContactCell(obj, vecMCell, lineDFM)
            iiCell = [];
            pf = lineDFM(1,:);
            for j=1:length(vecMCell)
                ID_MCell = vecMCell(j);
                vecPoint = obj.vecVert{ID_MCell};
                for i=1:length(vecPoint)
                    if norm(vecPoint(i,:)-pf) < 1e-10
                        iiCell = [iiCell, ID_MCell];
                    end
                end
            end
        end
        
        function [center, centerDFM, lineDFM, vertDFM, volFrac] ...
                =  calculateCenter2(obj, frac, iCellFrac1, iCellFrac2)
            ifrac1 = frac.vecIDFrac(iCellFrac1);
            ifrac2 = frac.vecIDFrac(iCellFrac2);
            vecIfrac = [ifrac1, ifrac1, ifrac2, ifrac2];
            Line1 = frac.vecCenterLineFrac{iCellFrac1};
            Line2 = frac.vecCenterLineFrac{iCellFrac2};          
            vertFrac = [Line1;Line2];
            % calculate the center of the intersection
            dx21 = vertFrac(2,1)-vertFrac(1,1);
            dx43 = vertFrac(4,1)-vertFrac(3,1);
            dy21 = vertFrac(2,2)-vertFrac(1,2);
            dy43 = vertFrac(4,2)-vertFrac(3,2);
            center = [0,0,vertFrac(1,3)];
            if abs(dx21) < 1e-10
                center(1:2) = [vertFrac(1,1), dy43/dx43*(vertFrac(1,1)-vertFrac(3,1))+vertFrac(3,2)];
            elseif abs(dx43) < 1e-10
                center(1:2) = [vertFrac(3,1), dy21/dx21*(vertFrac(3,1)-vertFrac(1,1))+vertFrac(1,2)];              
            elseif abs(dy21) < 1e-10
                center(1:2) = [dx43/dy43*(vertFrac(1,2)-vertFrac(3,2))+vertFrac(3,1), vertFrac(1,2)];
            elseif abs(dy43) < 1e-10
                center(1:2) = [dx21/dy21*(vertFrac(3,2)-vertFrac(1,2))+vertFrac(1,1), vertFrac(3,2)];
            else
                center(1) = (dx43*dy21*vertFrac(1,1)-dx21*dy43*vertFrac(3,1))/(dx43*dy21-dx21*dy43);
                center(2) = dy21/dx21*(center(1)-vertFrac(1,1))+vertFrac(1,2);
            end
                     
            volFrac = 0;
            centerDFM = (vertFrac+center)/2;   
            lineDFM = cell(0); % line passing through the center horizontally to the x-y plane.
            vertDFM = cell(0);
            dz = obj.dz(1, 1);            
            for i=1:4
                ifrac = vecIfrac(i);
                lineDFM_temp = [vertFrac(i,:); center];
                lineDFM(i) = {lineDFM_temp};
                
                if abs(vertFrac(i,2)-center(2)) < 1e-10                                
                    halfAper = frac.aperture(ifrac)/2; %/cos(obj.strike(ifrac)/180*pi);               
                    vertDFM_temp(1:2, :) = ones(2,1)*lineDFM_temp(1, :) + ...
                                                       [0,-halfAper,0;
                                                        0,halfAper,0;];                
                    vertDFM_temp(3:4, :) = ones(2,1)*lineDFM_temp(2, :) + ...
                                                       [0,-halfAper,0;
                                                        0,halfAper,0;];          
                else                           
                    halfAper = frac.aperture(ifrac)/2; %/cos(obj.strike(ifrac)/180*pi);                
                    vertDFM_temp(1:2, :) = ones(2,1)*lineDFM_temp(1, :) + ...
                                                       [-halfAper,0,0;
                                                        halfAper,0,0;];                 
                    vertDFM_temp(3:4, :) = ones(2,1)*lineDFM_temp(2, :) + ...
                                                       [-halfAper,0,0;
                                                        halfAper,0,0;];                 
                end 
                
                vertDFM_temp = sortrows(vertDFM_temp);
                tempVert = vertDFM_temp(4,:);
                vertDFM_temp(4,:) = vertDFM_temp(3,:);
                vertDFM_temp(3,:) = tempVert;
                vertDFM(i) = {vertDFM_temp}; % ld - lu - ru - rd
                
                volFrac(i) = dz*frac.aperture(ifrac)*norm(lineDFM_temp(1, :)-lineDFM_temp(2, :));
            end
                
            
        end
        
        %% Split a matrix into two parts
        function [obj, frac, connection] = splitMatrix(obj, frac, connection, ID_MCell_org, iCellFrac, ifrac)
            ID_MCell = ID_MCell_org;
            ID_MCell_1 = ID_MCell;
            ID_MCell_2 = ID_MCell + 1;         
            [c1, c2, elemVert1, elemVert2, vol1, vol2]...
                = obj.findSplitCenter(frac, ID_MCell, iCellFrac, ifrac);
            
            % Insert one element for each array
            [obj, frac] = obj.insertElement(frac, ID_MCell_1, ID_MCell_2, iCellFrac);            
            
            % Check boundary
            isBound = obj.isBoundary(1,ID_MCell_org);            
            min_x1 = min(elemVert1(:,1)); min_x2 = min(elemVert2(:,1));        
            max_x1 = max(elemVert1(:,1)); max_x2 = max(elemVert2(:,1));  
            vVert_org = obj.vecVert{ID_MCell_org};
            min_x = max(vVert_org(:,1)); max_x = max(vVert_org(:,1));
            obj.isBoundary(1,ID_MCell_1) = 0;
            obj.isBoundary(1,ID_MCell_2) = 0;
            if isBound == 1
                if min_x == min_x1
                    obj.isBoundary(1,ID_MCell_1) = 1;
                end
                if min_x == min_x2
                    obj.isBoundary(1,ID_MCell_2) = 1;  
                end        
            elseif isBound == 2
                if max_x == max_x1
                    obj.isBoundary(1,ID_MCell_1) = 2;
                end
                if max_x == max_x2
                    obj.isBoundary(1,ID_MCell_2) = 2;  
                end                  
            end
            
            % Update the matrix 1
            obj.volCell(1,ID_MCell_1) = vol1;
            obj.xCell(1,ID_MCell_1) = c1(1);
            obj.yCell(1,ID_MCell_1) = c1(2);
            obj.zCell(1,ID_MCell_1) = c1(3);
            obj.vecVert(ID_MCell_1) = {elemVert1};
            obj.nEmbeddedFrac(1,ID_MCell_1) = obj.nEmbeddedFrac(1,ID_MCell_1) + 1;            
            
            % Update the matrix 2
            obj.volCell(1,ID_MCell_2) = vol2;
            obj.xCell(1,ID_MCell_2) = c2(1);
            obj.yCell(1,ID_MCell_2) = c2(2);
            obj.zCell(1,ID_MCell_2) = c2(3);
            obj.vecVert(ID_MCell_2) = {elemVert2};
            obj.kx(1,ID_MCell_2) = obj.kx(1,ID_MCell_1);
            obj.ky(1,ID_MCell_2) = obj.ky(1,ID_MCell_1);
            obj.kz(1,ID_MCell_2) = obj.kz(1,ID_MCell_1);
            obj.poro(1,ID_MCell_2) = obj.poro(1,ID_MCell_1);
 
            obj.vGridIDToMultiFrac(ID_MCell_2) = {[obj.vGridIDToMultiFrac{ID_MCell_1}]};
            obj.vGridIDToMultiFCell(ID_MCell_2) = {[obj.vGridIDToMultiFCell{ID_MCell_1}]};
            obj.nEmbeddedFrac(1,ID_MCell_2) = obj.nEmbeddedFrac(1,ID_MCell_1);
            obj.vFracType(ID_MCell_2) = {[obj.vFracType{ID_MCell_1}]};
            obj.shapeFactor(1,ID_MCell_2) = obj.shapeFactor(1,ID_MCell_1);
            
            % Split fractures 
            [frac, obj] = frac.splitFracture(obj, ID_MCell_1, ID_MCell_2, iCellFrac);
            
            % Insert a new connection
            connection = connection.splitConnection(obj, frac, iCellFrac, ID_MCell_1, ID_MCell_2);
        end
        
        %% Split a matrix into two parts
        function [c1, c2, elemVert1, elemVert2, vol1, vol2] = findSplitCenter(obj, frac, ID_MCell, iCellFrac, ifrac)            
            c_m = [obj.xCell(1,ID_MCell), obj.yCell(1,ID_MCell), obj.zCell(1,ID_MCell)];
            lineEDFM = frac.vecCenterLineFrac{iCellFrac};  
            lineEDFM = calcIntersectionPoint(obj, obj.vecVert{ID_MCell}, lineEDFM(1,:), lineEDFM(2,:));
            cellVertex = [obj.vecVert{ID_MCell};lineEDFM];

            % sort vertEDFM (-dx/2,-dy/2,0)->(-dx/2,dy/2,0)->(dx/2,dy/2,0)
            % -> (dx/2,-dy/2,0) -> line 1 -> line 2
            dx = obj.dx(1, ID_MCell);
            dy = obj.dy(1, ID_MCell);
            dz = obj.dz(1, ID_MCell);

            % Chect the shape of split: triangle - penta or quad-quad
            isAbove = ((cellVertex-frac.center(:,ifrac)')*frac.normals(:,ifrac)>-1e-10);
            if abs(sum(isAbove(1:4))) == 2
                isTriangle = false;
            elseif abs(sum(isAbove(1:4))) == 1
                isTriangle = true;
            else
                isTriangle = true;
                zeroindex = (isAbove==1);
                isAbove(isAbove==0) = 1;
                isAbove(zeroindex) = 0;
            end
            isAbove(5:6) = 1;
            
            % Select vertices of one selected part
            elemVert1 = cellVertex(isAbove, :);            
            elemVert1 = obj.rearrangePoints(elemVert1);
            areaElem1 = 0;
            c_elem1 = 0;
            % Calculate a center mass of one part
            if isTriangle
                L1 = elemVert1(1,:)-elemVert1(2,:);
                L2 = elemVert1(1,:)-elemVert1(3,:);
                areaElem1 = norm(cross(L1, L2))/2;
                c_elem1 = mean(elemVert1,1);
            else
                avgPoint = mean(elemVert1,1);
                for i=1:4
                    c_elem(i,:) = elemVert1(i,:)+avgPoint;
                    L1 = elemVert1(i,:)-avgPoint;
                    if i == 4
                        L2 = elemVert1(1,:)-avgPoint;  
                        c_elem(i,:) = c_elem(i,:) + elemVert1(1,:);
                    else
                        L2 = elemVert1(i+1,:)-avgPoint;
                        c_elem(i,:) = c_elem(i,:) + elemVert1(i+1,:);
                    end
                    areaElem(i) = norm(cross(L1, L2))/2;
                    areaElem1 = areaElem1 + areaElem(i);
                    c_elem(i,:) = c_elem(i,:)/3;
                    c_elem1 = c_elem1 + areaElem(i) * c_elem(i,:);
                end    
            end            
            c_elem1 = c_elem1 / areaElem1;
            vol1 = areaElem1 * dz; 
            
            % Calculate the mass center of another part.
            zeroindex = (isAbove==1);
            isAbove(isAbove==0) = 1;
            isAbove(zeroindex) = 0;
            isAbove(5:6) = 1;
            elemVert2 = cellVertex(isAbove, :);
                        
            vVert = obj.vecVert{ID_MCell};           
            L1 = vVert(2,:)-vVert(1,:);             
            L2 = vVert(3,:)-vVert(1,:);  
            if length(vVert) == 4                
                L3 = vVert(4,:)-vVert(1,:);
                areaTot = 0.5*(abs(norm(cross(L1, L2)))...
                                +abs(norm(cross(L2, L3))));
            else
                areaTot = 0.5*abs(cross(L1, L2));                
            end
            areaElem2 = areaTot - areaElem1;
            c_elem2 = (c_m-areaElem1/areaTot*c_elem1)/(areaElem2/areaTot);
            vol2 = areaElem2 * dz; 
            
            c1 = c_elem1;
            c2 = c_elem2;
            elemVert1 = obj.rearrangePoints(elemVert1);
            elemVert2 = obj.rearrangePoints(elemVert2);
        end
        
        function elemVert = rearrangePoints(obj, elemVert)            
            elemVert = sortrows(elemVert); % ld - lu - rd - ru
            rowtemp = elemVert(3,:);
            elemVert(3,:) = elemVert(4,:);
            elemVert(4,:) = rowtemp; % ld - lu - ru - rd
        end
        
        function [obj, frac] = insertElement(obj, frac, ID_MCell_1, ID_MCell_2, iCellFrac)            
            iPrev = 1:ID_MCell_1; iNext = ID_MCell_2:length(obj.derivPoro);
            obj.volCell = [obj.volCell(1,iPrev), 0, obj.volCell(1,iNext)];
            obj.xCell = [obj.xCell(1,iPrev), 0, obj.xCell(1,iNext)];
            obj.yCell = [obj.yCell(1,iPrev), 0, obj.yCell(1,iNext)];
            obj.zCell = [obj.zCell(1,iPrev), 0, obj.zCell(1,iNext)];
            obj.vecVert = [obj.vecVert(1,iPrev), {[]}, obj.vecVert(1,iNext)];
            obj.dx = [obj.dx(1,iPrev), obj.dx(1,ID_MCell_1), obj.dx(1,iNext)];
            obj.dy = [obj.dy(1,iPrev), obj.dy(1,ID_MCell_1), obj.dy(1,iNext)];
            obj.dz = [obj.dz(1,iPrev), obj.dz(1,ID_MCell_1), obj.dz(1,iNext)];
            obj.kx = [obj.kx(1,iPrev), 0, obj.kx(1,iNext)];
            obj.ky = [obj.ky(1,iPrev), 0, obj.ky(1,iNext)];
            obj.kz = [obj.kz(1,iPrev), 0, obj.kz(1,iNext)];
            obj.poro = [obj.poro(1,iPrev), 0, obj.poro(1,iNext)];
            obj.poroPrev = [obj.poroPrev(1,iPrev), 0, obj.poroPrev(1,iNext)];
            obj.derivPoro = [obj.derivPoro(1,iPrev), 0, obj.derivPoro(1,iNext)];
            obj.nEmbeddedFrac = [obj.nEmbeddedFrac(1,iPrev), 0, obj.nEmbeddedFrac(1,iNext)];
            obj.vFracType = [obj.vFracType(1,iPrev), {[]}, obj.vFracType(1,iNext)];
            obj.shapeFactor = [obj.shapeFactor(1,iPrev), 0, obj.shapeFactor(1,iNext)];            
            obj.isBoundary =  [obj.isBoundary(1,iPrev), 0, obj.isBoundary(1,iNext)];
            obj.nMatrix = obj.nMatrix + 1;
            obj.nCell = obj.nCell + 1;
                           
            iCellFracIncrease = ID_MCell_2; 
            for i=1:length(obj.vGridIDToMultiFCell)
                vGridToFCell = [obj.vGridIDToMultiFCell{i}];
                increaseIndex = (vGridToFCell>=iCellFracIncrease);
                vGridToFCell(increaseIndex) = vGridToFCell(increaseIndex) + 1;
                obj.vGridIDToMultiFCell(i) = {vGridToFCell};
            end           
            
            for i=1:length(frac.vecIDCell)
                vecIDCell = frac.vecIDCell{i};
                increaseIndex = (vecIDCell>=iCellFracIncrease);
                vecIDCell(increaseIndex) = vecIDCell(increaseIndex) + 1;
                splitIndex = (vecIDCell == ID_MCell_1);
                if sum(splitIndex)==1
                    vecIDCell= [vecIDCell, (ID_MCell_1+1):ID_MCell_2];
                end
                frac.vecIDCell(i) = {unique(vecIDCell)};
            end                 
            
            obj.vGridIDToMultiFrac = [obj.vGridIDToMultiFrac(1,iPrev), {[]}, obj.vGridIDToMultiFrac(1,iNext)];
            obj.vGridIDToMultiFCell = [obj.vGridIDToMultiFCell(1,iPrev), {[]}, obj.vGridIDToMultiFCell(1,iNext)];            
            obj.mapOldToNew(obj.mapOldToNew>=ID_MCell_2) = obj.mapOldToNew(obj.mapOldToNew>=ID_MCell_2)+1;
        end
        
        function [obj, frac] = insertFracture(obj, frac, ID_MCell_1, ID_MCell_2, iCellFrac)            
            iPrev = 1:ID_MCell_1; iNext = ID_MCell_2:length(obj.derivPoro);
            obj.volCell = [obj.volCell(1,iPrev), 0, obj.volCell(1,iNext)];
            obj.xCell = [obj.xCell(1,iPrev), 0, obj.xCell(1,iNext)];
            obj.yCell = [obj.yCell(1,iPrev), 0, obj.yCell(1,iNext)];
            obj.zCell = [obj.zCell(1,iPrev), 0, obj.zCell(1,iNext)];
            obj.vecVert = [obj.vecVert(1,iPrev), {[]}, obj.vecVert(1,iNext)];
            obj.dx = [obj.dx(1,iPrev), obj.dx(1,ID_MCell_1), obj.dx(1,iNext)];
            obj.dy = [obj.dy(1,iPrev), obj.dy(1,ID_MCell_1), obj.dy(1,iNext)];
            obj.dz = [obj.dz(1,iPrev), obj.dz(1,ID_MCell_1), obj.dz(1,iNext)];
            obj.kx = [obj.kx(1,iPrev), 0, obj.kx(1,iNext)];
            obj.ky = [obj.ky(1,iPrev), 0, obj.ky(1,iNext)];
            obj.kz = [obj.kz(1,iPrev), 0, obj.kz(1,iNext)];
            obj.poro = [obj.poro(1,iPrev), 0, obj.poro(1,iNext)];
            obj.poroPrev = [obj.poroPrev(1,iPrev), 0, obj.poroPrev(1,iNext)];
            obj.derivPoro = [obj.derivPoro(1,iPrev), 0, obj.derivPoro(1,iNext)];
            obj.nEmbeddedFrac = [obj.nEmbeddedFrac(1,iPrev), 0, obj.nEmbeddedFrac(1,iNext)];
            obj.vFracType = [obj.vFracType(1,iPrev), {[]}, obj.vFracType(1,iNext)];
            obj.shapeFactor = [obj.shapeFactor(1,iPrev), 0, obj.shapeFactor(1,iNext)];            
            obj.isBoundary =  [obj.isBoundary(1,iPrev), 0, obj.isBoundary(1,iNext)];
            obj.nFrac = obj.nFrac + 1;
            obj.nCell = obj.nCell + 1;
            
            obj.vGridIDToMultiFrac = [obj.vGridIDToMultiFrac(1,iPrev), {[]}, obj.vGridIDToMultiFrac(1,iNext)];
            obj.vGridIDToMultiFCell = [obj.vGridIDToMultiFCell(1,iPrev), {[]}, obj.vGridIDToMultiFCell(1,iNext)];  
        end
        
        %% Define the first DPDP cell
        function obj = addDPDK(obj, frac, ID_MCell, ID_FCell, ifrac)            
            iCellFrac = ID_FCell - obj.nMatrix;
            obj.nCell = obj.nCell + 1; % add a fracture cell
            obj.nFrac = obj.nFrac + 1; % add a fracture cell 
            obj.nEmbeddedFrac(1,ID_MCell) = obj.nEmbeddedFrac(1,ID_MCell) + 1;
            obj.volCell(1, ID_FCell) = frac.vecVolFrac(1,iCellFrac);
            obj.kx(1, ID_FCell) = frac.fracPerm(1,ifrac);
            obj.ky(1, ID_FCell) = frac.fracPerm(1,ifrac);
            obj.kz(1, ID_FCell) = frac.fracPerm(1,ifrac);
            obj.poro(1, ID_FCell) = 1;
            obj.poroPrev(1, ID_FCell) = 0;
            obj.derivPoro(1, ID_FCell) = 0;            
            obj.xCell(1,ID_FCell) = obj.xCell(1,ID_MCell);
            obj.yCell(1,ID_FCell) = obj.yCell(1,ID_MCell);
            obj.zCell(1,ID_FCell) = obj.zCell(1,ID_MCell);  
            obj.vecVert(ID_FCell) = {[]};
            obj.shapeFactor(1, ID_FCell) = obj.calculateShapeFactor(obj.volCell(1, ID_FCell), ID_MCell);
            obj.isBoundary(ID_FCell) = 0;
            obj.nEmbeddedFrac(ID_FCell) = 0;
            obj.vFracType(ID_FCell) = {[]};
            obj.vGridIDToMultiFrac(ID_FCell) = {[]};
            obj.vGridIDToMultiFCell(ID_FCell) = {[]};
            obj.isBoundary(ID_FCell) = 0;            
        end
        
        %% Merge additional DPDP cells to the first DPDK cell.
        function [obj, frac] = addDPDKtoDPDK(obj, frac, ID_MCell, ID_FCell, ifrac, ID_FCell_rep)              
            vo11 = obj.volCell(1, ID_FCell_rep);
            vol2 = frac.vecVolFrac(1,ifrac);
            perm1 = obj.kx(1, ID_FCell_rep);
            perm2 = frac.fracPerm(1,ifrac);
            sh1 = obj.shapeFactor(1, ID_FCell_rep);
            sh2 = obj.calculateShapeFactor(vol2, ID_MCell);
            volTot = vo11+vol2;
            obj.kx(1, ID_FCell_rep) = vo11/volTot*perm1+vol2/volTot*perm2;
            obj.ky(1, ID_FCell_rep) = vo11/volTot*perm1+vol2/volTot*perm2;
            obj.kz(1, ID_FCell_rep) = vo11/volTot*perm1+vol2/volTot*perm2;
            obj.shapeFactor(1, ID_FCell_rep) = vo11/volTot*sh1+vol2/volTot*sh2;
            obj.volCell(1, ID_FCell_rep) = volTot;
            [frac, obj] = frac.removeFracCell(obj, ID_FCell); 
        end
        
        %% Calculate shape factor
        function shapeF = calculateShapeFactor(obj, volFrac, ID_MCell)
            dx = obj.dx(ID_MCell);
            dy = obj.dy(ID_MCell);
            dz = obj.dz(ID_MCell);
            p = [-1, dx+dy+dz, -(dx*dy+dy*dz+dz*dz), volFrac];
            Lfrac = roots(p);
            Lfrac = Lfrac(Lfrac>0);
            Lfrac = Lfrac(Lfrac < min([dx,dy,dz]));
            Lfrac = Lfrac(1);
            Lmx = dx - Lfrac;
            Lmy = dy - Lfrac;
            Lmz = dz - Lfrac;
            shapeF = 4*(1/Lmx^2) ; % Kazemi 
            %shapeF = 4*(1/Lmx^2+1/Lmy^2+1/Lmz^2) ; % Kazemi            
        end
        
        %% Calculate distance between EDFM and matrix
        function avgDist = calculateDistance(obj, iCellFrac, frac, ID_MCell)  
            ID_FCell = iCellFrac + obj.nMatrix;
            ifrac = frac.vecIDFrac(1,iCellFrac);
            
            c_m = [obj.xCell(1,ID_MCell), obj.yCell(1,ID_MCell), obj.zCell(1,ID_MCell)];
            lineEDFM = frac.vecCenterLineFrac{ID_FCell - obj.nMatrix};
            
            % sort vertEDFM (-dx/2,-dy/2,0)->(-dx/2,dy/2,0)->(dx/2,dy/2,0)
            % -> (dx/2,-dy/2,0) -> line 1 -> line 2
            dx = obj.dx(1, ID_MCell);
            dy = obj.dy(1, ID_MCell);
            dz = obj.dz(1, ID_MCell);
            cellVertex =[-dx,-dy,0;-dx,dy,0;dx,dy,0;dx,-dy,0]/2+c_m;
            cellVertex = [cellVertex;lineEDFM];
            
            % Chect the shape of split: triangle - penta or quad-quad
            isAbove = ((cellVertex-frac.center(:,ifrac)')*frac.normals(:,ifrac)>-1e-10);            
            isAbove(5:6) = 0;
            if abs(sum(isAbove)) == 2
                isTriangle = false;
            elseif abs(sum(isAbove)) == 1
                isTriangle = true;
            else
                isTriangle = true;
                zeroindex = (isAbove==1);
                isAbove(isAbove==0) = 1;
                isAbove(zeroindex) = 0;
            end
            isAbove(5:6) = 1;
            
            % Select vertices of one selected part
            elemVert = cellVertex(isAbove, :);
            elemVert = obj.rearrangePoints(elemVert);
            areaElem1 = 0;
            c_elem1 = 0;
            % Calculate a center mass of one part
            if isTriangle
                L1 = elemVert(1,:)-elemVert(2,:);
                L2 = elemVert(1,:)-elemVert(3,:);
                areaElem1 = norm(cross(L1, L2))/2;
                c_elem1 = mean(elemVert,1);
            else
                avgPoint = mean(elemVert,1);
                for i=1:4
                    c_elem(i,:) = elemVert(i,:)+avgPoint;
                    L1 = elemVert(i,:)-avgPoint;
                    if i == 4
                        L2 = elemVert(1,:)-avgPoint;  
                        c_elem(i,:) = c_elem(i,:) + elemVert(1,:);
                    else
                        L2 = elemVert(i+1,:)-avgPoint;
                        c_elem(i,:) = c_elem(i,:) + elemVert(i+1,:);
                    end
                    areaElem(i) = norm(cross(L1, L2))/2;
                    areaElem1 = areaElem1 + areaElem(i);
                    c_elem(i,:) = c_elem(i,:)/3;
                    c_elem1 = c_elem1 + areaElem(i) * c_elem(i,:);
                end    
            end            
            c_elem1 = c_elem1 / areaElem1;
            dist_1 = abs((c_elem1-frac.center(:,ifrac)')*frac.normals(:,ifrac))...
                     /norm(frac.normals(:,ifrac));
                 
            % Calculate the mass center of another part.
            areaTot = dx*dy;
            areaElem2 = areaTot - areaElem1;
            c_elem2 = (c_m-areaElem1/areaTot*c_elem1)/(areaElem2/areaTot);
            dist_2 = abs((c_elem2-frac.center(:,ifrac)')*frac.normals(:,ifrac))...
                         /norm(frac.normals(:,ifrac));
            avgDist = (dist_1*areaElem1+dist_2*areaElem2)/areaTot;
        end
        
        %% calculate an intersection point
        function vecInt = calcIntersectionPoint(obj, elemVert, f1, f2)
            nInt = 0;
            vecInt = zeros(0,3);
            for i=1:length(elemVert)
                % check intersection
                if i<length(elemVert)
                    basisEdge = elemVert(i+1,:)-elemVert(i,:);
                else
                    basisEdge = elemVert(1,:)-elemVert(i,:);               
                end
                basisEdge = [basisEdge(2), -basisEdge(1), 0];
                basisEdge = basisEdge/norm(basisEdge);
                innerP1 = (f1-elemVert(i,:))*basisEdge';
                innerP2 = (f2-elemVert(i,:))*basisEdge';
                
                if innerP1*innerP2 < 0
                    basisFrac = f1-f2;
                    basisFrac = basisFrac/norm(basisFrac);
                    basisFracEdge = (elemVert(i,:)-f2)/norm((elemVert(i,:)-f2));
                    distFrac2ToEdge = abs((elemVert(i,:)-f2)*basisEdge');
                    distFrac1Comp = abs((f1-f2)*basisEdge');
                    distChange = sign(basisFracEdge*basisFrac')...
                                  *distFrac2ToEdge/distFrac1Comp*(f1-f2);                    
                    nInt = nInt + 1;
                    vecInt(nInt,:) = f2+distChange;
                elseif abs(innerP1) < 1e-10
                    nInt = nInt + 1;
                    vecInt(nInt,:) = f1;
                elseif abs(innerP2) < 1e-10
                    nInt = nInt + 1;
                    vecInt(nInt,:) = f2;                    
                end
            end
        end
        
        %% Calculate distance between fractures
        function [area_avg, dist_avg, perm_avg, area_avg2, dist_avg2, perm_avg2] = ...
        calcIntersection2(obj, frac, ID_MCell, ID_FCell, ID_EDFM1, ID_EDFM2)
            ifrac = frac.vecIDFrac(1,ID_FCell-obj.nMatrix);
            ifrac_EDFM1 = frac.vecIDFrac(1,ID_EDFM1-obj.nMatrix); 
            ifrac_EDFM2 = frac.vecIDFrac(1,ID_EDFM2-obj.nMatrix);  
            lineEDFM = frac.vecCenterLineFrac{ID_FCell - obj.nMatrix};
            lineEDFM1 = frac.vecCenterLineFrac{ID_EDFM1 - obj.nMatrix};  
            lineEDFM2 = frac.vecCenterLineFrac{ID_EDFM2 - obj.nMatrix};          
            vertFrac = [lineEDFM;lineEDFM1;lineEDFM2];
            
            % calculate the center of the intersection
            dx21 = vertFrac(2,1)-vertFrac(1,1);
            dx43 = vertFrac(4,1)-vertFrac(3,1);
            dy21 = vertFrac(2,2)-vertFrac(1,2);
            dy43 = vertFrac(4,2)-vertFrac(3,2);
            center = [0,0,vertFrac(1,3)];
            if abs(dx21) < 1e-10
                center(1:2) = [vertFrac(1,1), dy43/dx43*(vertFrac(1,1)-vertFrac(3,1))+vertFrac(3,2)];
            elseif abs(dx43) < 1e-10
                center(1:2) = [vertFrac(3,1), dy21/dx21*(vertFrac(3,1)-vertFrac(1,1))+vertFrac(1,2)];              
            elseif abs(dy21) < 1e-10
                center(1:2) = [dx43/dy43*(vertFrac(1,2)-vertFrac(3,2))+vertFrac(3,1), vertFrac(1,2)];
            elseif abs(dy43) < 1e-10
                center(1:2) = [dx21/dy21*(vertFrac(3,2)-vertFrac(1,2))+vertFrac(1,1), vertFrac(3,2)];
            else
                center(1) = (dx43*dy21*vertFrac(1,1)-dx21*dy43*vertFrac(3,1))/(dx43*dy21-dx21*dy43);
                center(2) = dy21/dx21*(center(1)-vertFrac(1,1))+vertFrac(1,2);
            end
            
            dz = obj.dz(ID_MCell);
            dist_elem = zeros(6,1);
            T_elem = zeros(6,1);
            k_elem = [obj.kx(ID_FCell), obj.kx(ID_FCell),...
                      obj.kx(ID_EDFM1), obj.kx(ID_EDFM1),...
                      obj.kx(ID_EDFM2), obj.kx(ID_EDFM2)];
            A_elem = [frac.aperture(ifrac), frac.aperture(ifrac), ...
                      frac.aperture(ifrac_EDFM1), frac.aperture(ifrac_EDFM1), ...
                      frac.aperture(ifrac_EDFM2), frac.aperture(ifrac_EDFM2)]*dz;     
            for i=1:6
                % Ti = ki*Ai/Di
                dist_elem(i) = norm(center-vertFrac(i,:))/2;    
                if dist_elem(i) == 0
                    T_elem(i) = 0;
                else
                    T_elem(i) = k_elem(i)*A_elem(i)/dist_elem(i);
                end
            end            
            T_sum = sum(T_elem);
            %T_12 = (T1T3+T1T4+T2T3+T2T4)/(T1+T2+T3+T4)
            T_12 = sum(T_elem(1:2))*sum(T_elem(3:4))/T_sum;
            area_avg = obj.dz(ID_MCell)*frac.aperture(ifrac);
            dist_avg = dist_elem(1)+dist_elem(2);
            % T_12 = perm_avg * A_avg / dist_avg
            perm_avg = T_12/area_avg*dist_avg; 
            
            %T_12 = (T1T3+T1T4+T2T3+T2T4)/(T1+T2+T3+T4)
            T_34 = sum(T_elem(3:4))*sum(T_elem(5:6))/T_sum;
            area_avg2 = obj.dz(ID_MCell)*frac.aperture(ifrac_EDFM1);
            dist_avg2 = dist_elem(3)+dist_elem(4);
            % T_12 = perm_avg * A_avg / dist_avg
            perm_avg2 = T_34/area_avg2*dist_avg2; 
            
        end
        
        %% Calculate distance between fractures
        function [area_avg, dist_avg, perm_avg] = calculateIntersectionDistance(obj, frac, ID_MCell, ID_FCell, ID_FCell_embed)
            ifrac = frac.vecIDFrac(1,ID_FCell-obj.nMatrix);
            ifrac_embed = frac.vecIDFrac(1,ID_FCell_embed-obj.nMatrix);  
            lineEDFM = frac.vecCenterLineFrac{ID_FCell - obj.nMatrix};
            lineEDFM_embed = frac.vecCenterLineFrac{ID_FCell_embed - obj.nMatrix};            
            vertFrac = [lineEDFM;lineEDFM_embed];
            
            % calculate the center of the intersection
            dx21 = vertFrac(2,1)-vertFrac(1,1);
            dx43 = vertFrac(4,1)-vertFrac(3,1);
            dy21 = vertFrac(2,2)-vertFrac(1,2);
            dy43 = vertFrac(4,2)-vertFrac(3,2);
            center = [0,0,vertFrac(1,3)];
            if abs(dx21) < 1e-10
                center(1:2) = [vertFrac(1,1), dy43/dx43*(vertFrac(1,1)-vertFrac(3,1))+vertFrac(3,2)];
            elseif abs(dx43) < 1e-10
                center(1:2) = [vertFrac(3,1), dy21/dx21*(vertFrac(3,1)-vertFrac(1,1))+vertFrac(1,2)];              
            elseif abs(dy21) < 1e-10
                center(1:2) = [dx43/dy43*(vertFrac(1,2)-vertFrac(3,2))+vertFrac(3,1), vertFrac(1,2)];
            elseif abs(dy43) < 1e-10
                center(1:2) = [dx21/dy21*(vertFrac(3,2)-vertFrac(1,2))+vertFrac(1,1), vertFrac(3,2)];
            else
                center(1) = (dx43*dy21*vertFrac(1,1)-dx21*dy43*vertFrac(3,1))/(dx43*dy21-dx21*dy43);
                center(2) = dy21/dx21*(center(1)-vertFrac(1,1))+vertFrac(1,2);
            end
            
            dz = obj.dz(ID_MCell);
            dist_elem = zeros(4,1);
            T_elem = zeros(4,1);
            k_elem = [obj.kx(ID_FCell), obj.kx(ID_FCell),...
                      obj.kx(ID_FCell_embed), obj.kx(ID_FCell_embed)];
            A_elem = [frac.aperture(ifrac), frac.aperture(ifrac), ...
                      frac.aperture(ifrac_embed), frac.aperture(ifrac_embed)]*dz;
            for i=1:4
                % Ti = ki*Ai/Di
                dist_elem(i) = norm(center-vertFrac(i,:))/2;    
                if dist_elem(i) == 0
                    T_elem(i) = 0;
                else
                    T_elem(i) = k_elem(i)*A_elem(i)/dist_elem(i);
                end
            end            
            T_sum = sum(T_elem);
            %T_12 = (T1T3+T1T4+T2T3+T2T4)/(T1+T2+T3+T4)
            T_12 = (T_elem(1)+T_elem(2))*(T_elem(3)+T_elem(4))/T_sum;
            area_avg = obj.dz(ID_MCell)*frac.aperture(ifrac);
            dist_avg = dist_elem(1)+dist_elem(2);
            % T_12 = perm_avg * A_avg / dist_avg
            perm_avg = T_12/area_avg*dist_avg; 
        end
        
        %% Calculate distance between fractures
        function [area_avg, dist_avg, perm_avg] = calculateIntersectionDFM(obj, frac, ID_MCell, ID_FCell, ID_FCell_embed)
            ifrac = frac.vecIDFrac(1,ID_FCell-obj.nMatrix);
            ifrac_embed = frac.vecIDFrac(1,ID_FCell_embed-obj.nMatrix);  
            lineEDFM = frac.vecCenterLineFrac{ID_FCell - obj.nMatrix};
            lineEDFM_embed = frac.vecCenterLineFrac{ID_FCell_embed - obj.nMatrix}; 
            iCellFrac = ID_FCell-obj.nMatrix;
            
            center =  frac.vecCenterIntersect{iCellFrac};
            ID_F_inter = frac.vecListIntersect{iCellFrac}+obj.nMatrix;
            iCellFrac_inter = frac.vecListIntersect{iCellFrac};
            ifrac_inter = frac.vecIDFrac(iCellFrac_inter);
            vertFrac = zeros(0,3);
            for i=1:length(ID_F_inter)
                vertFrac = [vertFrac;frac.vecCenterLineFrac{iCellFrac_inter(i)}];
            end
            dz = obj.dz(ID_MCell);
            dist_elem = zeros(4,1);
            T_elem = zeros(4,1);
            k_elem = obj.kx(ID_F_inter);
            A_elem = frac.aperture(ifrac_inter)*dz;
            for i=1:4
                % Ti = ki*Ai/Di
                dist_elem(i) = norm(center-vertFrac(2*i-1,:))/2;    
                if dist_elem(i) == 0
                    T_elem(i) = 0;
                else
                    T_elem(i) = k_elem(i)*A_elem(i)/dist_elem(i);
                end
                
                if ID_FCell == ID_F_inter(i)
                    T1 = T_elem(i);
                    D1 = dist_elem(i);
                elseif ID_FCell_embed == ID_F_inter(i)
                    T2 = T_elem(i);
                    D2 = dist_elem(i);
                end
            end            
            T_sum = sum(T_elem);
            %T_12 = (T1T3+T1T4+T2T3+T2T4)/(T1+T2+T3+T4)
            T_12 = T1*T2/T_sum;
            area_avg = obj.dz(ID_MCell)*frac.aperture(ifrac);
            dist_avg = D1+D2;
            % T_12 = perm_avg * A_avg / dist_avg
            perm_avg = T_12/area_avg*dist_avg; 
        end
        
        %% Update porosity with respect to Po
        function arg = updatePoro(obj, Poil)
            arg = obj.poro_ref*exp(obj.compRock*(Poil-obj.poro_P_ref));
        end
        
        %% Updated(Porosity)/d(Po) at (i,j,k) at time n+1
        function arg = updateDerivPoro(obj, Poil)
            arg = obj.poro_ref*obj.compRock*exp(obj.compRock*(Poil-obj.poro_P_ref));
        end
        
        %% Backup porosity
        function [poroPrev] = backUpGrid(obj)            
            poroPrev = obj.poro;
        end
        
        %% Implement all update functions in the class
        function [arg1, arg2] = updateGrid(obj, fluid)            
            arg1 = updatePoro(obj, fluid.Po);
            arg2 = updateDerivPoro(obj, fluid.Po);
        end
    end
end

