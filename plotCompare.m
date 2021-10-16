function [iFigure] = plotCompare(iFigure, outputRef, output, fluidRef, fluid, wellRef, well, grid, fn)
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
plot(outputRef.simTime/365, outputRef.WOPR);
hold on
plot(output.simTime/365, output.WOPR);
hold on
plot(time/365, WOPR, '--g');
hold on
lgd = legend('Reference', 'New well location', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('Well oil rate(STB/day)','FontSize',14)
title('Well oil rate at the production well','FontSize',14)
set(gca,'FontSize',14)

%% 2. Plot BHP at production well
subplot(2,2,2);
plot(outputRef.simTime/365, outputRef.WBHPPROD);
hold on
plot(output.simTime/365, output.WBHPPROD);
hold on
plot(time/365, WBHPPROD, '--g');
hold on
lgd = legend('Reference', 'New well location', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('Production well BHP(psi)','FontSize',14)
title('Bottomhole pressure at the production well','FontSize',14)
set(gca,'FontSize',14)
iFigure = iFigure + 1;

%% 3. Plot Water rate at production well
subplot(2,2,3);
plot(outputRef.simTime/365, outputRef.WWPR);
hold on
plot(output.simTime/365, output.WWPR);
hold on
plot(time/365, WWPR, '--g');
hold on
lgd = legend('Reference', 'New well location', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('Water rate at the production well[STB/day]','FontSize',14)
title('Water rate at the production well','FontSize',14)
set(gca,'FontSize',14)

%% 4.BHP at injection well
subplot(2,2,4);
plot(outputRef.simTime/365, outputRef.WBHPINJ);
hold on
plot(output.simTime/365, output.WBHPINJ);
hold on
plot(time/365, WBHPINJ, '--g');
hold on
lgd = legend('Reference', 'New well location', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('BHP at the injection well[PSI]','FontSize',14)
title('BHP at the injection well','FontSize',14)
set(gca,'FontSize',14)
iFigure = iFigure + 1;

%% 5. Plot the saturation field.
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
for i = 1:length(output.simTime)
    if output.simTime(i)<=(iT-1)*output.simTime(end)/4
        plotTime(iT) = output.simTime(i);
    else
        iT = iT + 1;
        plotTime(iT) = output.simTime(i);
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

%% 6. Plot saturation
for i = 1:size(output.SWAT,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(output.SWAT(i,:,:),grid.nz,grid.nx));
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
ss = suptitle('Water saturation with different well location');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

for i = 1:size(outputRef.SWAT,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(outputRef.SWAT(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinjRef, yinjRef,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprodRef, yprodRef,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime2(i)/365); 
    title(titleName, 'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Water saturation with original well location');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% 7. Plot porosity
for i = 1:size(output.PORO,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(output.PORO(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime(i)/365); 
    title(titleName, 'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Porosity with compressibility');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

for i = 1:size(outputRef.PORO,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(outputRef.PORO(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinjRef, yinjRef,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprodRef, yprodRef,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime2(i)/365); 
    title(titleName, 'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Porosity with incompressibility');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% 8. Plot pressure
for i = 1:size(output.PRESSURE,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(output.PRESSURE(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinj, yinj,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprod, yprod,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime(i)/365); 
    title(titleName, 'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Oil pressure with different well location');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

for i = 1:size(outputRef.PRESSURE,1)
    figure(iFigure)
    subplot(2,3,i);
    imagesc(reshape(outputRef.PRESSURE(i,:,:),grid.nz,grid.nx));
    colormap jet
    colorbar
    hold on
    plot(xinjRef, yinjRef,'marker','o','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    plot(xprodRef, yprodRef,'marker','x','MarkerSize',8, 'Color', 'r', 'LineWidth',2)
    titleName = sprintf('Year %.1f', plotTime2(i)/365); 
    title(titleName, 'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Oil pressure with original well location');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% 9. Plot cummulative oil production
figure(iFigure)
plot(outputRef.simTime/365, outputRef.FOPT);
hold on
plot(output.simTime/365, output.FOPT);
lgd = legend('Reference', 'New well location');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('Cummulative oil production(STB)','FontSize',14)
title('Cummulative oil production','FontSize',14)
set(gca,'FontSize',14)
end
