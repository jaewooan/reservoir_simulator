function [iFigure] = plotCompare3wells(iFigure, outputRef, output, fluidRef, fluid, wellRef, well, grid, fn)
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
    time(index) = str2double(splitWord{1,iSpace+1}); %%days
    WOPR(index) = str2double(splitWord{1,iSpace+3}); %% well oil rate(STB/DAY) 
    WWPR(index) = str2double(splitWord{1,iSpace+4}); %% well water rate(STB/DAY)
    WWIR(index) = str2double(splitWord{1,iSpace+5}); %% Water injection rate(Psia)
    WBHPPROD(index) = str2double(splitWord{1,iSpace+6}); %% Prod BHP(Psia)
    WBHPINJ(index) = str2double(splitWord{1,iSpace+7}); %% Inj BHP(Psia)
    FPR(index) = str2double(splitWord{1,iSpace+8}); %% Average reservoir pressure(Psia)
    FOPT(index) = str2double(splitWord{1,iSpace+9}); %% Cumulative oil production of the field(STB)
    FWPT(index) = str2double(splitWord{1,iSpace+10}); %% cumulative water production of the field(STB)
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
    FOE(index) = str2double(splitWord{1,iSpace+2}); %% request oil recovery 
    FOPV(index) = str2double(splitWord{1,iSpace+3}); %%field oil pore volume(RB)   
    FWPV(index) = str2double(splitWord{1,iSpace+4}); %%field water pore volume(rb)
    BPR_prod(index) = str2double(splitWord{1,iSpace+5}); %%block pressure requested(psia)
    BPR_inj(index) = str2double(splitWord{1,iSpace+6}); %%block pressure requested(psia)
    BPR_prod2(index) = str2double(splitWord{1,iSpace+7}); %%block pressure requested(psia)
    WOPR2(index) = str2double(splitWord{1,iSpace+8}); %% well oil rate(STB/DAY) 
    WWPR2(index) = str2double(splitWord{1,iSpace+9}); %% well water rate(STB/DAY)
    WBHPPROD2(index) = str2double(splitWord{1,iSpace+10}); %% Prod BHP(Psia)
    tline = fgetl(fid);
    index = index + 1;
end

totWOPR = 0; totWWPR = 0;
for i = 1:well.numberOfWells
    if well.wellType(i)==1
        totWOPR = totWOPR+output.WOPR(i,:);
        totWWPR = totWWPR+output.WWPR(i,:);
    end
end

%% 1. Plot oil rate
figure(iFigure)
subplot(2,2,1);
plot(outputRef.simTime/365, outputRef.WOPR, 'r');
hold on
plot(output.simTime/365, output.WOPR(1,:), 'k', 'LineWidth', 2);
hold on
plot(output.simTime/365, output.WOPR(2,:), 'g', 'LineWidth', 2);
hold on
plot(output.simTime/365, totWOPR, 'b');
hold on
plot(time/365, WOPR, 'm--');
hold on
plot(time/365, WOPR2, 'm--');
hold on
lgd = legend('Case 3', 'Prod 1', 'Prod 2','Prod 1+2', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('Well oil rate(STB/day)','FontSize',14)
title('Well oil rate at the production well','FontSize',14)
set(gca,'FontSize',14)

%% 2. Plot BHP at production well
subplot(2,2,2);
plot(outputRef.simTime/365, outputRef.WBHPPROD, 'r');
hold on
plot(output.simTime/365, output.WBHPPROD(1,:),'b');
hold on
plot(output.simTime/365, output.WBHPPROD(2,:),'k');
hold on
plot(time/365, WBHPPROD, 'm--');
hold on
plot(time/365, WBHPPROD2, 'm--');
hold on
lgd = legend('Case 3', 'Prod 1', 'Prod 2', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('Production well BHP(psi)','FontSize',14)
title('Bottomhole pressure at the production well','FontSize',14)
set(gca,'FontSize',14)
iFigure = iFigure + 1;

%% 3. Plot Water rate at production well
subplot(2,2,3);
plot(outputRef.simTime/365, outputRef.WWPR, 'r');
hold on
plot(output.simTime/365, output.WWPR(1,:), 'k', 'LineWidth', 2);
hold on
plot(output.simTime/365, output.WWPR(2,:), 'g', 'LineWidth', 2);
hold on
plot(output.simTime/365, totWWPR, 'b');
hold on
plot(time/365, WWPR, 'm--');
hold on
plot(time/365, WWPR2, 'm--');
hold on
lgd = legend('Case 3', 'Prod 1', 'Prod 2', 'Prod 1+2', 'Eclipse');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('Water rate at the production well[STB/day]','FontSize',14)
title('Water rate at the production well','FontSize',14)
set(gca,'FontSize',14)

%% 4.BHP at injection well
subplot(2,2,4);
plot(outputRef.simTime/365, outputRef.WBHPINJ,'b');
hold on
plot(output.simTime/365, output.WBHPINJ,'k');
hold on
plot(time/365, WBHPINJ, 'm--');
hold on
lgd = legend('Case 3', 'Case 5', 'Eclpise');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('BHP at the injection well[PSI]','FontSize',14)
title('BHP at the injection well','FontSize',14)
set(gca,'FontSize',14)
iFigure = iFigure + 1;

%% Load the location of wells
iProdCellRef = wellRef.wellRevMap(1);
iInjCellRef = wellRef.wellRevMap(2);
xprodRef = mod((iProdCellRef-1),30)+1;
yprodRef = floor((iProdCellRef-1)/30)+1;
xinjRef = mod((iInjCellRef-1),30)+1;
yinjRef = floor((iInjCellRef-1)/30)+1;

iProdCell = well.wellRevMap(1);
iProd2Cell = well.wellRevMap(2);
iInjCell = well.wellRevMap(3);
xprod = mod((iProdCell-1),30)+1;
yprod = floor((iProdCell-1)/30)+1;
xprod2 = mod((iProd2Cell-1),30)+1;
yprod2 = floor((iProd2Cell-1)/30)+1;
xinj = mod((iInjCell-1),30)+1;
yinj = floor((iInjCell-1)/30)+1;


%% Plot saturation
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
    imagesc(reshape(output.SWAT(i,:,:),grid.nz,grid.nx));
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
    plot(xprod2, yprod2,'marker','x','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
end
xlabel('x(cell #)');
ylabel('z(cell #)');
set(h,'FontSize',14);    
colormap jet
hold on;
ss = suptitle('Water saturation');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% Plot porosity
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
    hold on;
    plot(xprod2, yprod2,'marker','x','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
    titleName = sprintf('Year %.1f', plotTime(i)/365); 
    title(titleName, 'FontSize',14);
    set(gca,'FontSize',14)
    hold on
end
ss = suptitle('Water saturation with different well loaction');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% Plot pressure
figure(iFigure)
for i = 1:5%size(output.PRESSURE,1)  
    sp = subplot(2, 3, i);
    imagesc(reshape(output.PRESSURE(i,:,:),grid.nz,grid.nx));
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
    plot(xprod2, yprod2,'marker','x','MarkerSize',8, 'Color', 'w', 'LineWidth',2)
    hold on;
end
ss = suptitle('Oil pressure');
set(ss,'FontSize',20,'FontWeight','normal')
iFigure=iFigure+1;

%% Plot cummulative oil production
figure(iFigure)
plot(outputRef.simTime/365, outputRef.FOPT);
hold on
plot(output.simTime/365, output.FOPT);
lgd = legend('Case 3', 'Case 5');
lgd.FontSize = 14;
xlabel('Time(year)','FontSize',14)
ylabel('Cummulative oil production(STB)','FontSize',14)
title('Cummulative oil production','FontSize',14)
set(gca,'FontSize',14)
end
