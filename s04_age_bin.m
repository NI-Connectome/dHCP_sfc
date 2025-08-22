% 3.3 Overgrowth of structure-function coupling in preterm infants 
% bin infants into several groups and generate figure 3c 

clc; clear; close all

load ./data/NetsInfo.mat
col = [0.627450980392157,0,0.0705882352941177;
    1,0.431372549019608,0.294117647058824;
    0.454901960784314,0.827450980392157,1];

for dataset = 3

switch dataset
    case 3
        load ./data/dHCP/FT_436.mat
    case 1
        load ./data/dHCP/PT1_120.mat
    case 2
        load ./data/dHCP/PT2_117.mat
end
pma = info.scan_age;
whR = mean(R);

wdow_range = 20;                                          
wdow_nums = 10;                                           
step = ceil((length(pma)-wdow_range)/wdow_nums);          
wdow_tail = wdow_range:step:length(pma);                  

%% bins
[sortedAge, idx] = sort(pma);
clear cage age_bins_whR age_bins_ntR

for w = 1:wdow_nums
    tailw = wdow_tail(w);
    idx_s = idx(tailw-wdow_range+1:tailw); 
    cage(w) = mean(pma(idx_s));
    age_bins_whR(w) = mean(whR(idx_s));
    age_bins_std(w) = std(whR(idx_s));
end

% figure('Position',[50 50 600 450]);
ysm = smooth(cage',age_bins_whR',8,'rlowess'); 
age_bins_std = smooth(cage',age_bins_std',8,'rlowess'); 
hold on 

%% Add standard deviation intervals
xconf = [cage, cage(end:-1:1)];
yconf = [ysm'-age_bins_std', ysm(end:-1:1)'+age_bins_std'];
p = fill(xconf,yconf, col(dataset,:));
p.FaceAlpha = 0.5;
p.EdgeColor = 'none';
uistack(p, 'down')

Line = plot(cage,ysm,'color',col(dataset,:));
% errorbar(cage, age_bins_wBrain, error(1,:), 'k.', 'LineWidth', 1); 
set(Line, 'LineStyle', '-','LineWidth', 2);
% set(Line, 'Marker', 'd','MarkerEdgeColor',[0,0,0],'MarkerFaceColor','w','MarkerSize',6);


end

ylabel('Whole brain coupling index');
xlabel('PMA(weeks)');

legend('Preterm at TEA','Full-term','Location', 'north','Orientation','horizon')
legend boxoff
set(gca,'TickDir', 'out', 'TickLength', [.02 .02],'FontWeight', 'bold',...
    'fontname','arial','LineWidth',1.7, 'Fontsize',16);
xlim([min(cage)-0.2, max(cage)+0.2])