function [iFigure] = plotCompare3(iFigure, outputRef, outputSmall, outputLarge,...
    fluidRef, fluidSmall, fluidLarge, wellRef, well, grid, fn1, fn2, fn3)
%% Read a reference file
fid = fopen(strcat(pwd, '/reference/',fn1,'.RSM'), 'r');
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
    tline = fgetl(fid);
    index = index + 1;
end

fid2 = fopen(strcat(pwd, '/reference/',fn2,'.RSM'), 'r');
tline = fgetl(fid2);
nCount = 0;
while nCount < 3
    if contains(tline, '--')
        nCount = nCount + 1;
    end    
    tline = fgetl(fid2);
end
index = 1;
while ~contains(tline, '1                 ')
    splitWord = strsplit(tline);
    iSpace = isempty(splitWord{1,1});
    timeS(index) = str2double(splitWord{1,iSpace+1}); %days
    WOPRS(index) = str2double(splitWord{1,iSpace+3}); % well oil rate(STB/DAY) 
    WWPRS(index) = str2double(splitWord{1,iSpace+4}); % well water rate(STB/DAY)
    WWIRS(index) = str2double(splitWord{1,iSpace+5}); % Water injection rate(Psia)
    WBHPPRODS(index) = str2double(splitWord{1,iSpace+6}); % Prod BHP(Psia)
    WBHPINJS(index) = str2double(splitWord{1,iSpace+7}); % Inj BHP(Psia)
    tline = fgetl(fid2);
    index = index + 1;
end

fid3 = fopen(strcat(pwd, '/reference/',fn3,'.RSM'), 'r');
tline = fgetl(fid3);
nCount = 0;
while nCount < 3
    if contains(tline, '--')
        nCount = nCount + 1;
    end    
    tline = fgetl(fid3);
end
index = 1;
while ~contains(tline, '1                 ')
    splitWord = strsplit(tline);
    iSpace = isempty(splitWord{1,1});
    timeL(index) = str2double(splitWord{1,iSpace+1}); %days
    WOPRL(index) = str2double(splitWord{1,iSpace+3}); % well oil rate(STB/DAY) 
    WWPRL(index) = str2double(splitWord{1,iSpace+4}); % well water rate(STB/DAY)
    WWIRL(index) = str2double(splitWord{1,iSpace+5}); % Water injection rate(Psia)
    WBHPPRODL(index) = str2double(splitWord{1,iSpace+6}); % Prod BHP(Psia)
    WBHPINJL(index) = str2double(splitWord{1,iSpace+7}); % Inj BHP(Psia)
    tline = fgetl(fid3);
    index = index + 1;
end


%% 1. Plot oil rate
figure(iFigure)
subplot(2,2,1);
plot(outputSmall.simTime/365, outputSmall.WOPR,'r', 'LineWidth', 2);
hold on
plot(outputRef.simTime/365, outputRef.WOPR,'b', 'LineWidth', 2);
hold on
plot(outputLarge.simTime/365, outputLarge.WOPR,'k', 'LineWidth', 2);
hold on
plot(time/365, WOPR, '--g');
hold on
plot(timeS/365, WOPRS, '--g');
hold on
plot(timeL/365, WOPRL, '--g');
hold on
lgd = legend('Small', 'Base', 'Large', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(day)','FontSize',14)
ylabel('Well oil rate(STB/day)','FontSize',14)
title('Well oil rate at the production well','FontSize',14)
set(gca,'FontSize',14)

%% 2. Plot Water rate at production well
subplot(2,2,2);
plot(outputSmall.simTime/365, outputSmall.WWPR,'r', 'LineWidth', 2);
hold on
plot(outputRef.simTime/365, outputRef.WWPR,'b', 'LineWidth', 2);
hold on
plot(outputLarge.simTime/365, outputLarge.WWPR,'k', 'LineWidth', 2);
hold on
plot(time/365, WWPR, '--g');
hold on
plot(timeS/365, WWPRS, '--g');
hold on
plot(timeL/365, WWPRL, '--g');
hold on
lgd = legend('Small', 'Base', 'Large', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(day)','FontSize',14)
ylabel('Water rate[STB/day]','FontSize',14)
title('Water rate at the production well','FontSize',14)
set(gca,'FontSize',14)

%% 3.BHP at injection well
subplot(2,2,3);
plot(outputSmall.simTime/365, outputSmall.WBHPINJ,'r', 'LineWidth', 2);
hold on
plot(outputRef.simTime/365, outputRef.WBHPINJ,'b', 'LineWidth', 2);
hold on
plot(outputLarge.simTime/365, outputLarge.WBHPINJ,'k', 'LineWidth', 2);
hold on
plot(time/365, WBHPINJ, '--g');
hold on
plot(timeS/365, WBHPINJS, '--g');
hold on
plot(timeL/365, WBHPINJL, '--g');
hold on
lgd = legend('Small', 'Base', 'Large', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(day)','FontSize',14)
ylabel('BHP at the injection well[PSI]','FontSize',14)
title('BHP at the injection well','FontSize',14)
set(gca,'FontSize',14)
iFigure = iFigure + 1;

%% 5. Check well location
iProdCellRef = wellRef.wellRevMap(1);
iInjCellRef = wellRef.wellRevMap(2);
xprodRef = mod((iProdCellRef-1),30)+1;
yprodRef = floor((iProdCellRef-1)/30)+1;
xinjRef = mod((iInjCellRef-1),30)+1;
yinjRef = floor((iInjCellRef-1)/30)+1;

iProdCell = well.wellRevMap(1);
iInjCell = well.wellRevMap(2);
xprod = mod((iProdCell-1),30)+1;
yprod = floor((iProdCell-1)/30)+1;
xinj = mod((iInjCell-1),30)+1;
yinj = floor((iInjCell-1)/30)+1;

plotTime = zeros(1,5); iT = 2;
for i = 1:length(outputSmall.simTime)
    if outputSmall.simTime(i)<=(iT-1)*outputSmall.simTime(end)/4
        plotTime(iT) = outputSmall.simTime(i);
    else
        iT = iT + 1;
        plotTime(iT) = outputSmall.simTime(i);
    end    
end
plotTime = floor(plotTime/1);

plotTime2 = zeros(1,5); iT = 2;
for i = 1:length(outputRef.simTime)
    if outputRef.simTime(i)<=(iT-1)*outputRef.simTime(end)/4
        plotTime2(iT) = outputRef.simTime(i);
    else
        iT = iT + 1;
        plotTime2(iT) = outputRef.simTime(i);
    end    
end
plotTime2 = floor(plotTime2/1);

plotTime3 = zeros(1,5); iT = 2;
for i = 1:length(outputLarge.simTime)
    if outputLarge.simTime(i)<=(iT-1)*outputLarge.simTime(end)/4
        plotTime3(iT) = outputLarge.simTime(i);
    else
        iT = iT + 1;
        plotTime3(iT) = outputLarge.simTime(i);
    end    
end
plotTime3 = floor(plotTime3/1);

%% 6. Plot saturation
for i = 1:size(outputSmall.SWAT,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(outputSmall.SWAT(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime(i)/365); 
    title(titleName,'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Water saturation with small mobility');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;


for i = 1:size(outputRef.SWAT,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(outputRef.SWAT(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime2(i)/365); 
    title(titleName,'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Water saturation with base mobility');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

for i = 1:size(outputLarge.SWAT,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(outputLarge.SWAT(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime3(i)/365); 
    title(titleName,'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Water saturation with large mobility');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% 6. Plot saturation
for i = 1:size(outputSmall.PRESSURE,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(outputSmall.PRESSURE(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime(i)/365); 
    title(titleName,'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Oil pressure with small mobility');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;


for i = 1:size(outputRef.PRESSURE,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(outputRef.PRESSURE(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime2(i)/365); 
    title(titleName,'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Oil pressure with base mobility');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

for i = 1:size(outputLarge.PRESSURE,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(outputLarge.PRESSURE(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime3(i)/365); 
    title(titleName,'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Oil pressure with large mobility');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% Plot cummulative oil production
figure(iFigure)
plot(outputSmall.simTime/365, outputSmall.FOPT,'r');
hold on
plot(outputRef.simTime/365, outputRef.FOPT,'b');
hold on
plot(outputLarge.simTime/365, outputLarge.FOPT,'k');
lgd = legend('Small', 'Base', 'Large');
lgd.FontSize = 14;
xlabel('Time(day)','FontSize',14)
ylabel('Cummulative oil production(STB)','FontSize',14)
title('Cummulative oil production','FontSize',14)
set(gca,'FontSize',14)
end
