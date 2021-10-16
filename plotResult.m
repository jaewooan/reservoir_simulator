function [iFigure] = plotResult(iFigure, output, fn, fluid, grids, well, frac,...
                                output2, fn2, fluid2, grids2, well2, frac2)
%% Read a reference file
fid = fopen(strcat(pwd, '/reference/',fn,'.RSM'), 'r');

tline = fgetl(fid);
nCount = 0;
while nCount < 3
    if contains(tline, '--')
        nCount = nCount + 1;
    end    
    tline = fgetl(fid);
end

index = 1;
while ~contains(tline, '1                 ')
    splitWord = strsplit(tline);
    iSpace = isempty(splitWord{1,1});
    time(index) = str2double(splitWord{1,iSpace+1}); %days
    WOPR(index) = str2double(splitWord{1,iSpace+3}); % well oil rate(STB/DAY) 
    WWPR(index) = str2double(splitWord{1,iSpace+4}); % well water rate(STB/DAY)
    WWIR(index) = str2double(splitWord{1,iSpace+5}); % Water injection rate(Psia)
    WBHPPROD(index) = str2double(splitWord{1,iSpace+6}); % Prod BHP(Psia)
    WBHPINJ(index) = str2double(splitWord{1,iSpace+7}); % Inj BHP(Psia)
    FPR(index) = str2double(splitWord{1,iSpace+8}); % Average reservoir pressure(Psia)
    FOPT(index) = str2double(splitWord{1,iSpace+9}); % Cumulative oil production of the field(STB)
    FWPT(index) = str2double(splitWord{1,iSpace+10}); % cumulative water production of the field(STB)
    tline = fgetl(fid);
    index = index + 1;
end

nCount = 0;
while nCount < 3
    if contains(tline, '--')
        nCount = nCount + 1;
    end    
    tline = fgetl(fid);
end

index = 1;
while tline > -1
    splitWord = strsplit(tline);
    iSpace = isempty(splitWord{1,1});
    FOE(index) = str2double(splitWord{1,iSpace+2}); % request oil recovery 
    FOPV(index) = str2double(splitWord{1,iSpace+3}); %field oil pore volume(RB)   
    FWPV(index) = str2double(splitWord{1,iSpace+4}); %field water pore volume(rb)
    BPR_prod(index) = str2double(splitWord{1,iSpace+5}); %block pressure requested(psia)
    BPR_inj(index) = str2double(splitWord{1,iSpace+6}); %block pressure requested(psia)
    tline = fgetl(fid);
    index = index + 1;
end

%% 1. Plot oil rate
figure(iFigure)
subplot(2,2,1);
%plot(time/365, WOPR);
%hold on
plot(output.simTime/365, output.WOPR);
lgd = legend('Eclipse', 'Simulator');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('Well oil rate(STB/day)','FontSize',14)
title('Well oil rate at the production well','FontSize',14)
set(gca,'FontSize',14)

%% 2. Plot BHP at production well
subplot(2,2,2);
%plot(time/365, WBHPPROD);
%hold on
%plot(output.simTime/365, output.WBHPPROD);
plot(output.simTime/365, output.WWPR);
lgd = legend('Eclipse', 'Simulator');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('BHP at the production well(psi)','FontSize',14)
title('BHP at the production well','FontSize',14)
set(gca,'FontSize',14)

%% 3.BHP at injection well
subplot(2,2,3);
plot(time/365, WBHPINJ);
hold on
plot(output.simTime/365, output.WBHPINJ);
lgd = legend('Eclipse', 'Simulator');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('BHP at the injection well[PSI]','FontSize',14)
title('BHP at the injection well','FontSize',14)
set(gca,'FontSize',14)
iFigure = iFigure + 1;

%% 4. Read RSM files
if contains(fn, "60yrs")
    refOutput = read_ecl(strcat(pwd, '/reference/',fn,'.X0010'));    
else
    refOutput = read_ecl(strcat(pwd, '/reference/',fn,'.X0010'));    
end
for row = 1:grids.ny
    for col = 1:grids.nx
%         if refOutput.SWAT(col+30*(row-1)) < 0
%             Sw_ref(row,col) = 0;
%         else
%             Sw_ref(row,col) = abs(refOutput.SWAT(col+30*(row-1)));
%         end
%         P_ref(row, col)= refOutput.PRESSURE(col+30*(row-1));
%         P_ref(row, col)= refOutput.PRESSURE(col+30*(row-1));
        % iCell = col+(row-1)*grids.nx;
        iCell = grids.mapPlotToGrid(col, row, 1);
%        permX(row,col) = grids.kx(iCell);
%         if grids.nEmbeddedFrac(iCell) > 0
%             IFCell = grids.vGridIDToMultiFCell{iCell};
%             permX(row,col) = grids.kx(IFCell);
%         end
    end
end

xprod = 24;%mod((iProdCell-1),grids.nx)+1;
yprod = 9;%floor((iProdCell-1)/grids.nx)+1;
xinj = 2;%mod((iInjCell-1),grids.nx)+1;
yinj = 15;%floor((iInjCell-1)/grids.nx)+1;



%% 5. Saturation of simulation
% Check plot time
plotTime = zeros(1,5); iT = 2;
for i = 1:length(output.simTime)
    if output.simTime(i)<=(iT-1)*output.simTime(end)/4
        plotTime(iT) = output.simTime(i);
    else
        iT = iT + 1;
        plotTime(iT) = output.simTime(i);
    end    
end
plotTime = floor(plotTime/1);

figure(iFigure)
for i = 1:5%size(output.SWAT,1)  
    sp = subplot(2, 3, i);
    ifraccell = (grids.ny*grids.nx+1):size(output.SWAT,2);
    for row = 1:grids.ny
        for col = 1:grids.nx
            iCell = grids.mapPermPlotToGrid(col, row, 1);
            SWATtemp(row,col) = output.SWAT(i,iCell);
        end
    end
%     SWATtemp(frac.vecIDCell(1,ifraccell(1:frac.vecNumCell(1))-grids.nMatrix))...
%         = output.SWAT(i, ifraccell(1:frac.vecNumCell(1)));
%   SWATtemp(frac.vecIDCell(1,ifraccell((frac.vecNumCell(1)+1):frac.nTotCell)-grids.nMatrix))...
%       = output.SWAT(i, ifraccell((frac.vecNumCell(1)+1):frac.nTotCell));
    %imagesc(reshape(output.SWAT(i,1:grids.ny*grids.nx), grids.ny, grids.nx)');
    %imagesc(reshape(SWATtemp, grids.ny, grids.nx)');
    imagesc(SWATtemp)
    set(gcf, 'Color', [1,1,1]); % white background
    h=colorbar;
    set(gca, 'FontSize', 14, 'YDir','reverse'); % bigger font size
    titleName = sprintf('Year %.1f', plotTime(i)/365); 
    title(titleName);
    xlabel('x(cell #)');
    ylabel('y(cell #)');
    set(h,'FontSize',14);    
    colormap jet
    hold on;
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
end
% subplot(2, 3, 6);
% imagesc(Sw_ref);
% set(gcf, 'Color', [1,1,1]); % white background
% h=colorbar;
% set(gca, 'FontSize', 14, 'YDir','reverse'); % bigger font size
% titleName = sprintf('Eclipse at year 5'); 
% title(titleName);
% xlabel('x(cell #)');
% ylabel('z(cell #)');
% set(h,'FontSize',14);    
% colormap jet
% hold on;
plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
hold on;
plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
hold on;
ss = suptitle('Water saturation');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% 5. Saturation of simulation
% Check plot time
plotTime2 = zeros(1,5); iT = 2;
for i = 1:length(output2.simTime)
    if output2.simTime(i)<=(iT-1)*output2.simTime(end)/4
        plotTime2(iT) = output2.simTime(i);
    else
        iT = iT + 1;
        plotTime2(iT) = output2.simTime(i);
    end    
end
plotTime2 = floor(plotTime2/1);

figure(iFigure)
for i = 1:5%size(output.SWAT,1)  
    sp = subplot(2, 3, i);
    ifraccell = (grids2.ny*grids2.nx+1):size(output2.SWAT,2);
    for row = 1:grids2.ny
        for col = 1:grids2.nx
            iCell = grids2.mapPermPlotToGrid(col, row, 1);
            SWATtemp2(row,col) = output2.SWAT(i,iCell);
        end
    end
%     SWATtemp(frac.vecIDCell(1,ifraccell(1:frac.vecNumCell(1))-grids.nMatrix))...
%         = output.SWAT(i, ifraccell(1:frac.vecNumCell(1)));
%   SWATtemp(frac.vecIDCell(1,ifraccell((frac.vecNumCell(1)+1):frac.nTotCell)-grids.nMatrix))...
%       = output.SWAT(i, ifraccell((frac.vecNumCell(1)+1):frac.nTotCell));
    %imagesc(reshape(output.SWAT(i,1:grids.ny*grids.nx), grids.ny, grids.nx)');
    %imagesc(reshape(SWATtemp, grids.ny, grids.nx)');
    imagesc(SWATtemp2)
    set(gcf, 'Color', [1,1,1]); % white background
    h=colorbar;
    set(gca, 'FontSize', 14, 'YDir','reverse'); % bigger font size
    titleName = sprintf('Year %.1f', plotTime2(i)/365); 
    title(titleName);
    xlabel('x(cell #)');
    ylabel('y(cell #)');
    set(h,'FontSize',14);    
    colormap jet
    hold on;
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
end
% subplot(2, 3, 6);
% imagesc(Sw_ref);
% set(gcf, 'Color', [1,1,1]); % white background
% h=colorbar;
% set(gca, 'FontSize', 14, 'YDir','reverse'); % bigger font size
% titleName = sprintf('Eclipse at year 5'); 
% title(titleName);
% xlabel('x(cell #)');
% ylabel('z(cell #)');
% set(h,'FontSize',14);    
% colormap jet
% hold on;
plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
hold on;
plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
hold on;
ss = suptitle('Water saturation');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% 6. Pressure DFM
figure(iFigure)
for i = 5:5%size(output.PRESSURE,1)  
   % sp = subplot(2, 3, i);
    
    for row = 1:grids.ny
        for col = 1:grids.nx
            iCell = grids.mapPermPlotToGrid(col, row, 1);
            Prestemp(row,col) = output.PRESSURE(i,iCell);
        end
    end
    imagesc(Prestemp);
    %imagesc(reshape(output.PRESSURE(i,1:grids.ny*grids.nx),grids.ny,grids.nx)');
    set(gcf, 'Color', [1,1,1]); % white background
    h=colorbar;
    set(gca, 'FontSize', 14, 'YDir','reverse'); % bigger font size
    ss = title('Oil Pressure with DFM (psi)');    
    set(ss,'FontSize',20,'FontWeight','normal')
    xlabel('x(cell #)');
    ylabel('y(cell #)');
    set(h,'FontSize',14);  
    
    r = [1 0 0];       %# start
    w = [.9 .9 .9];    %# middle
    b = [0 0 1];       %# end

    %# colormap of size 64-by-3, ranging from red -> white -> blue
    c1 = zeros(32,3); c2 = zeros(32,3);
    for i=1:3
        c1(:,i) = linspace(b(i), w(i), 32);
        c2(:,i) = linspace(w(i), r(i), 32);
    end
    c = [c1(1:end-1,:);c2];
    colormap(c)
    hold on;
end
iFigure=iFigure+1;

%% EDFM
figure(iFigure)
for i = 5:5%size(output.PRESSURE,1)      
    for row = 1:grids2.ny
        for col = 1:grids2.nx
            iCell = grids2.mapPermPlotToGrid(col, row, 1);
            Prestemp2(row,col) = output2.PRESSURE(i,iCell);
        end
    end
    imagesc(Prestemp2);
    set(gcf, 'Color', [1,1,1]); % white background
    h=colorbar;
    set(gca, 'FontSize', 14, 'YDir','reverse'); % bigger font size
    ss = title('Oil pressure with EDFM (psi)');
    set(ss,'FontSize',20,'FontWeight','normal')
    xlabel('x(cell #)');
    ylabel('y(cell #)');
    set(h,'FontSize',14);  
    
    r = [1 0 0];       %# start
    w = [.9 .9 .9];    %# middle
    b = [0 0 1];       %# end

    %# colormap of size 64-by-3, ranging from red -> white -> blue
    c1 = zeros(32,3); c2 = zeros(32,3);
    for i=1:3
        c1(:,i) = linspace(b(i), w(i), 32);
        c2(:,i) = linspace(w(i), r(i), 32);
    end
    c = [c1(1:end-1,:);c2];
    colormap(c)
    hold on;
end
iFigure=iFigure+1;

%% EDFM-DFM
figure(iFigure)
for i = 5:5%size(output.PRESSURE,1)  
    disp('err: ')
    sum(sum(abs(Prestemp2-Prestemp)/(grids.ny*grids.nx)))
    Pressdiff = (Prestemp2-Prestemp)./Prestemp*100;
    imagesc(Pressdiff);
    set(gcf, 'Color', [1,1,1]); % white background
    h=colorbar;
    set(gca, 'FontSize', 14, 'YDir','reverse'); % bigger font size
    ss = title('Relative error (%) in oil pressure');
    set(ss,'FontSize',20,'FontWeight','normal')
    xlabel('x(cell #)');
    ylabel('y(cell #)');
    set(h,'FontSize',14);  
    mean(mean(Pressdiff))
    r = [1 0 0];       %# start
    w = [.9 .9 .9];    %# middle
    b = [0 0 1];       %# end

    %# colormap of size 64-by-3, ranging from red -> white -> blue
    c1 = zeros(32,3); c2 = zeros(32,3);
    for i=1:3
        c1(:,i) = linspace(b(i), w(i), 32);
        c2(:,i) = linspace(w(i), r(i), 32);
    end
    c = [c1(1:end-1,:);c2];
    colormap(c)
    hold on;
end
iFigure=iFigure+1;

%% 8. Time step
figure(iFigure)
subplot(2,2,1);
plot(output.simTime/365, output.dT_out)
hold on
xlabel('Time(year)','FontSize',14);
ylabel('Time step(day)','FontSize',14);
set(gca,'FontSize',14);
title('Time step[day]','FontSize',14)

%% 9. Number of newton iteration
subplot(2,2,2);
plot(output.simTime/365, output.numIter)
hold on
xlabel('Time(year)','FontSize',14);
ylabel('# of iterations','FontSize',14);
set(gca,'FontSize',14);
title('The number of newton iteration','FontSize',14)

%% 10. Max CFL number
subplot(2,2,3);
plot(output.simTime/365, output.CFLmax)
hold on
xlabel('Time(year)','FontSize',14);
ylabel('Maximum CFL number','FontSize',14);
set(gca,'FontSize',14);
title('Maximum CFL number','FontSize',14);

%% 11. Convergence criteria
subplot(2,2,4);
plot(output.simTime/365, log10(output.maxR),'k')
hold on
plot(output.simTime/365, log10(output.maxS),'b')
hold on
plot(output.simTime/365, log10(output.maxP),'g')
hold on
plot([0,output.simTime(end)/365],[-3, -3], 'r--')
lgd = legend('Max of R', 'Sw change', 'Po change', 'Criteria');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14);
ylabel('Dimensionless value for each criteria','FontSize',14);
set(gca,'FontSize',14);
title('Convergence criteria','FontSize',14)
iFigure = iFigure+1;

%% 12. CFL number
figure(iFigure)
for i = 1:size(output.CFL,1)  
    sp = subplot(2, 3, i);
    imagesc(reshape(output.CFL(i,1:grids.ny*grids.nx),grids.ny,grids.nx)');
    set(gcf, 'Color', [1,1,1]); % white background
    h=colorbar;
    set(gca, 'FontSize', 14, 'YDir','reverse'); % bigger font size
    titleName = sprintf('Year %.1f', plotTime(i)/365); 
    title(titleName);
    xlabel('x(cell #)');
    ylabel('z(cell #)');
    set(h,'FontSize',14);    
    colormap jet
    hold on;
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
end
ss = suptitle('CFL number');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

end
