function [fluid_, grid_, connection_, equation_, well_, output_, frac_] = initialize(filename)
%% Read input files
filename = strcat(pwd, '/input/',filename);
[fluid_, grid_, connection_, equation_, well_, output_, frac_] =...
    readInputFile(filename); %p_initial [psi]
 
% Update pressure considering gravity effect.
%[fluid_.vecP, fluid_.Po]=fluid_.initializePressure(grid_);

% Initialize porosity; need to update for frac
[grid_.poro, grid_.derivPoro] = grid_.updateGrid(fluid_);

% Initialize fluid properties
[fluid_.bw, fluid_.bo, fluid_.Bw, fluid_.Bo, fluid_.gamma_denW, fluid_.gamma_denO,...
 fluid_.visW, fluid_.visO, fluid_.Krw, fluid_.Kro, fluid_.deriv_bw, fluid_.deriv_bo,...
 fluid_.derivVisW, fluid_.derivVisO, fluid_.derivKrw, fluid_.derivKro]...
 = fluid_.updateFluid();

% Initialize connection properties => update for M-F: always use matrix
% property for permeability.
connection_.transGeo = connection_.initializeTransGeo();
[connection_.iUpwind, connection_.transTot, connection_.avgDensity,...
    connection_.derivTransTot, connection_.derivAvgDensity]...
                = connection_.updateConnnection(fluid_, grid_);

% Initialize well transmissibility
[well_.kWB, well_.rWB, well_.WI] = well_.initializeWellTransGeo(grid_);
well_.Twell = well_.updateWellTrans(fluid_);

% Initialize variables at time level n and iteration k.
[fluid_.vecPPrev, fluid_.PoPrev, fluid_.SwPrev,...
    fluid_.bwPrev, fluid_.boPrev, fluid_.PoPrevIter, ...
    fluid_.SwPrevIter] = fluid_.backUpFluid();
[grid_.poroPrev] = grid_.backUpGrid();

% Grids for plot
grid_ = grid_.mappingGridtoPlot();

%% Plot perm
iFigure = 1;
xprod = 24;%mod((iProdCell-1),grid.nx)+1;
yprod = 9;%floor((iProdCell-1)/grid.nx)+1;
xinj = 2;%mod((iInjCell-1),grid.nx)+1;
yinj = 15;%floor((iInjCell-1)/grid.nx)+1;
for row = 1:grid_.ny
    for col = 1:grid_.nx
        %iCell = col+(row-1)*grid_.nx;
        iCell = grid_.mapPermPlotToGrid(col, row, 1);
        try
        permX(row,col) = grid_.kx(iCell);
        catch
            disp('dd')
        end
%         if grid_.nEmbeddedFrac(iCell) > 0
%             IFCell = grid_.vGridIDToMultiFCell{iCell};
%             permX(row,col) = grid_.kx(IFCell);
%         end
    end
end
figure(iFigure)
imagesc(permX);
colormap jet
colorbar
hold on
% plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
% hold on;
% plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
% hold on
grid on;
xlabel('x(cell #)','FontSize',14);
ylabel('y(cell #)','FontSize',14);
set(gca, 'FontSize', 14, 'YDir','reverse'); % bigger font size
title('Absolute permeability[mD]','FontSize',14)
end

