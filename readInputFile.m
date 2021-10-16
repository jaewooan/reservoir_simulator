function [fluid_, grid_, connection_, equation_, well_, output_, frac_] = ...    
        readInputFile(filename)

fid = fopen(filename, 'r' );
% Save grid information into a grid variable.
dd =1;
%% Grid
tline = fgetl(fid);
while ischar(tline) & (~contains(tline, "Fluid") && ~contains(tline, "Fractures"))  
    if (tline > -1) & contains(tline, "Nx")
        splitWord = strsplit(tline);
        nx = str2double(splitWord{1,3});
    end
        
    if (tline > -1) & contains(tline, "Ny")
        splitWord = strsplit(tline);
        ny = str2double(splitWord{1,3});
    end
    
    if (tline > -1) & contains(tline, "Nz")
        splitWord = strsplit(tline);
        nz = str2double(splitWord{1,3});
    end     
    
    if (tline > -1) & contains(tline, "Dx")
        splitWord = strsplit(tline);
        dx = str2double(splitWord{1,3})/dd;
    end
        
    if (tline > -1) & contains(tline, "Dy")
        splitWord = strsplit(tline);
        dy = str2double(splitWord{1,3})/dd;
    end
    
    if (tline > -1) & contains(tline, "Dz")
        splitWord = strsplit(tline);
        dz = str2double(splitWord{1,3})/dd;
    end    
      
    if (tline > -1) & contains(tline, "Poro")
        splitWord = strsplit(tline);
        poro = str2double(splitWord{1,3});
    end
        
    if (tline > -1) & contains(tline, "PermX")
        splitWord = strsplit(tline);
        kx = str2double(splitWord{1,3});
    end
    
    if (tline > -1) & contains(tline, "PermY")
        splitWord = strsplit(tline);
        ky = str2double(splitWord{1,3});
    end
        
    if (tline > -1) & contains(tline, "PermZ")
        splitWord = strsplit(tline);
        kz = str2double(splitWord{1,3});
    end
    
    if (tline > -1) & contains(tline, "Tops")
        splitWord = strsplit(tline);
        tops = str2double(splitWord{1,3});
    end
    
    if (tline > -1) & contains(tline, "RockComp")
        splitWord = strsplit(tline);
        rc = str2double(splitWord{1,3});
    end
    
    tline = fgetl(fid);
end

grid_ = grids(nx, ny, nz);
grid_.dx = dx*ones(1, grid_.nx*grid_.ny*grid_.nz);
grid_.dy = dy*ones(1, grid_.nx*grid_.ny*grid_.nz);
grid_.dz = dz*ones(1, grid_.nx*grid_.ny*grid_.nz);
grid_.poro_ref = poro;
grid_.kx = kx*ones(1, grid_.nx*grid_.ny*grid_.nz);
grid_.ky = ky*ones(1, grid_.nx*grid_.ny*grid_.nz);
grid_.kz = kz*ones(1, grid_.nx*grid_.ny*grid_.nz);
grid_.z_top = tops;
grid_.compRock = rc;

% Calculate basic properties of grid.
grid_ = grid_.initializeGrid();

% End grid.

%% Make a connection object
connection_=connection(grid_);
% End connection.

%% Read fracture information
isfrac = false;
frac_ = -1;
if contains(tline, "Fractures")
    isfrac = true;
    nfrac=0;
    fracCenter = zeros(3,1);
    fracLength = 0;
    fracWidth = 0;
    fracStrike = 0;
    fracAperture = 0;
    fracModel = 0;
    fracRangeX = 0;
    fracRangeY = 0;
    fracPoro = 0;
    fracPerm = 0;
    fracConductivity = 0;
    FLX=0;
    FLY=0;
    FLZ=0;
    while ischar(tline) && ~contains(tline, "Fluid")
        splitWord = strsplit(tline);
        if contains(tline, "label")    
            nfrac = str2double(splitWord{1,3});
        elseif contains(tline, "center") 
            fracCenter(1,nfrac) = str2double(splitWord{1,4})/dd;
            fracCenter(2,nfrac) = str2double(splitWord{1,5})/dd;
            fracCenter(3,nfrac) = dz*nz/2+tops;%str2double(splitWord{1,6});
        elseif contains(tline, "length")  
            fracLength(1,nfrac) = str2double(splitWord{1,3})/dd;
        elseif contains(tline, "width") 
            fracWidth(1,nfrac) = str2double(splitWord{1,3})/dd;
        elseif contains(tline, "rangeX") 
            fracRangeX(1,nfrac) = str2double(splitWord{1,3})/2/dd;
        elseif contains(tline, "rangeY") 
            fracRangeY(1,nfrac) = str2double(splitWord{1,3})/2/dd;
            fracRangeZ(1,nfrac) = dz*nz/2;
            fracStrike(1,nfrac) = 0;
            fracWidth(1,nfrac) = 0;
            fracLength(1,nfrac) = 0;
        elseif contains(tline, "strike") 
            fracStrike(1,nfrac) = str2double(splitWord{1,3});
            fracRangeX(1,nfrac) = 0;
            fracRangeY(1,nfrac) = 0;
            fracRangeZ(1,nfrac) = 0;
            FLX(1,nfrac) = 0;
            FLY(1,nfrac) = 0;
            FLZ(1,nfrac) = 0;
        elseif contains(tline, "aperture")
            fracAperture(1,nfrac) = str2double(splitWord{1,3})/1000000/0.3048; %micrometer->ft
            fracPoro(1,nfrac) = 0;        
            fracPerm(1,nfrac) = 0;
        elseif contains(tline, "permeability")
            fracPerm(1,nfrac) = str2double(splitWord{1,3});  
        elseif contains(tline, "model")
            if contains(tline, "DPDK")
                fracModel(1,nfrac)  = 0;            
            elseif contains(tline, "EDFM")
                fracModel(1,nfrac)  = 1;          
            elseif contains(tline, "DDFM")
                fracModel(1,nfrac)  = 2;          
            else
                fracModel(1,nfrac)  = 0;          
            end        
        end        
        tline = fgetl(fid);
    end
    % Define fracture
    frac_=fracture(fracCenter, fracLength, fracWidth, fracStrike,...
                   fracAperture, nfrac, fracModel, fracRangeX, fracRangeY,...
                   fracRangeZ, fracPoro, fracPerm, grid_);

    % Find location of cells with fractures
    [frac_, grid_] = frac_.initializeFracGeometry(grid_);

    % Embed fractures into grids
    [grid_, frac_, connection_] = grid_.defineFractureCells(frac_, connection_);
    
    % Update connection lists
    % at first, signle dpdk, signle edfm, intersection edfm, dpdk-edfm in one
    % cell, -> single dfm, intersection dfm, dpdk-dfm, dfm-edfm, dfm-edfm-dpdk
    connection_ = connection_.makeConnection(grid_, frac_);
end
% End fracture.

%% Create objects
equation_=equation(grid_, connection_);
fluid_=fluid(grid_);
output_=output(grid_);
well_=well(grid_);

%% Read fluid information
while ischar(tline) & ~contains(tline, "Wells")   
    tline = fgetl(fid);
    while ischar(tline) & ~contains(tline, "Wells") & ~contains(tline, "Water")
        splitWord = strsplit(tline);
        if contains(tline, "Comp")  
            fluid_.compOil = str2double(splitWord{1,3});    
        elseif contains(tline, "Vis")
            fluid_.visO_std = str2double(splitWord{1,3});
            fluid_.compVisOil = 0;
        elseif contains(tline, "Den")
            fluid_.denO_std = str2double(splitWord{1,3});            
        elseif contains(tline, "pres")
            fluid_.Po = str2double(splitWord{1,3})*ones(1, grid_.nCell);
            fluid_.vecP(1:2:end) = fluid_.Po;            
        end
        tline = fgetl(fid);
    end   
    
    while ischar(tline) & ~contains(tline, "Wells") & ~contains(tline, "Oil")
        splitWord = strsplit(tline);
        if contains(tline, "Comp")
            fluid_.compWat = str2double(splitWord{1,3});              
        elseif contains(tline, "Vis")
            fluid_.visW_std = str2double(splitWord{1,3}); 
            fluid_.compVisWat = 0;           
        elseif contains(tline, "Den")
            fluid_.denW_std = str2double(splitWord{1,3});            
        elseif contains(tline, "Sat")
            fluid_.Sw = str2double(splitWord{1,3})*ones(1, grid_.nCell);            
        end
        tline = fgetl(fid);
    end
end
% Finish to read fluid information

%% Read well information
nProd = 0; nInj = 0;
while ischar(tline) & ~contains(tline, "Boundary")   
    tline = fgetl(fid);
    while ischar(tline) & ~contains(tline, "well:") & ~contains(tline, "Boundary")     
        splitWord = strsplit(tline);
        if contains(tline, "type")            
            well_.numberOfWells = well_.numberOfWells + 1;  
            if contains(tline, "PROD")
                nProd = nProd+1;
                well_.wellType(well_.numberOfWells) = 1;
            else
                nInj = nInj+1;
                well_.wellType(well_.numberOfWells) = 2;
            end
        elseif contains(tline, "ii")   
            wellX = str2double(splitWord{1,3}); 
        elseif contains(tline, "jj")  
            wellY = str2double(splitWord{1,3});
        elseif contains(tline, "kk1")  
            wellZ = str2double(splitWord{1,3});  
            cellIndex = wellX + (wellY-1)*grid_.nx+(wellZ-1)*grid_.nx*grid_.ny;
            cellIndex = grid_.mapOldToNew(cellIndex);
            well_.wellMap(1,cellIndex) = well_.numberOfWells;
            well_.wellRevMap(1,well_.numberOfWells) = cellIndex;
            well_.Pwf(1,well_.numberOfWells) = 0;
            well_.qo(1,well_.numberOfWells) = 0;
            well_.qw(1,well_.numberOfWells) = 0;
        elseif contains(tline, "welldia")
            well_.wellDia(well_.numberOfWells) = str2double(splitWord{1,3});
            well_.wellR(well_.numberOfWells) = well_.wellDia(well_.numberOfWells)/2;
        elseif contains(tline, "rate")
            well_.flowConst(1,well_.numberOfWells) = str2double(splitWord{1,3});    
        elseif contains(tline, "bhp")
            well_.BHPconst(1,well_.numberOfWells) = str2double(splitWord{1,3});    
        elseif contains(tline, "mode:")
            if contains(tline, "RATE")
                well_.wellMode(1,well_.numberOfWells) = 1;
            else
                well_.wellMode(1,well_.numberOfWells) = 2;
            end     
        elseif contains(tline, "fluid")
            if contains(tline, "OIL")
                well_.qo(1,well_.numberOfWells) = well_.flowConst(well_.numberOfWells); 
            elseif contains(tline, "WAT")                
                well_.qw(1,well_.numberOfWells) = well_.flowConst(well_.numberOfWells); 
            end
        elseif contains(tline, "label")
            well_.wellName(well_.numberOfWells) = splitWord{1,3};           
        end
        tline = fgetl(fid); 
    end
end

while ischar(tline) & ~contains(tline, "Simulator") 
    tline = fgetl(fid);    
    splitWord = strsplit(tline);
    if (tline > -1) & contains(tline, "PL")
        connection_.PL = str2double(splitWord{1,3});
    elseif (tline > -1) & contains(tline, "PR")
        connection_.PR = str2double(splitWord{1,3}); 
    end
end

while ischar(tline)   
    tline = fgetl(fid);
    if (tline > -1) & contains(tline, "time")   
        splitWord = strsplit(tline);
        equation_.totT = 365*str2double(splitWord{1,3});
    end
end



fclose(fid);


end

