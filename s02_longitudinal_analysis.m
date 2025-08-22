% 3.2 Longitudinal analysis of structure-function 
% coupling for preterm infants
clc; clear; close all
addpath ./funct

%% data loading
load ./data/dHCP/FT_436.mat R
mFT = mean(mean(R));
load ./data/dHCP/PT1_120.mat info R
info1 = info;
PT1_origin = R;
load ./data/dHCP/PT2_117.mat info R
info2 = info; 
PT2_origin = R; 

%% find individuals with logitudinal data
cosub = intersect(info1.participant_id, info2.participant_id);
for s = 1:size(cosub)
    sub_idx(s,1) = find(info1.participant_id == cosub(s));
    sub_idx(s,2) = find(info2.participant_id == cosub(s));
end

%% data Extracting from two scans with the same subject
PT1 = PT1_origin(:,sub_idx(:,1));
PT2 = PT2_origin(:,sub_idx(:,2));
info_pt1 = info1(sub_idx(:,1),:);
info_pt2 = info2(sub_idx(:,2),:);

%% increased subjs
Ins = PT2-PT1;
Ins = Ins>0;
Ins_rate = sum(Ins,2)./63;

%% Paired t-tests

% whole brain
whR_PT1 = mean(PT1);
whR_PT2 = mean(PT2);
[h_wh,p_wh] = ttest(whR_PT1',whR_PT2');

% net
ntR_PT1 = netCI(PT1);
ntR_PT2 = netCI(PT2);

for n = 1:8
    x = ntR_PT1(n,:);
    y = ntR_PT2(n,:);
    [h_net(n,1),p_net(n,1)] = ttest(x',y');
end
fixedp_net = mafdr(p_net,'BHFDR',true);
fixedh_net = fixedp_net<0.05;

% ROI
for r = 1:90
    x = PT1(r,:);
    y = PT2(r,:);
    [h_roi(r,1),p_roi(r,1)] = ttest(x',y');
end
fixedp_roi = mafdr(p_roi,'BHFDR',true);
fixedh_roi = fixedp_roi<0.05;

%% individual analysis (Figure 2c)
mPT1 = mean(PT1);
mPT2 = mean(PT2);
     
[mPT1, idx] = sort(mPT1);

mPT2 = mPT2(idx);

figure
hold on
x = 1:length(mPT1);
for i = 1:length(x)
    if mPT1(i) - mPT2(i)<0
        col = addcolor(139);
    else col = addcolor(108);
    end
    plot([i,i], [mPT1(i), mPT2(i)], 'color', col,'linewidth',1)
end

a1 = scatter(x, mPT1,25, [160, 0, 18]./255,'filled')
a2 = scatter(x, mPT2,25, [255, 110, 75]./255,'filled')
a3 = plot([x(1)-2,x(end)+2], [mFT,mFT], '--', 'color', addcolor(200),...
    'Linewidth', 2)
hold off

box off
ax = gca;
set(gca,'TickDir', 'out', 'TickLength', [.02 .02],'FontWeight', 'bold',...
    'fontname','arial','LineWidth',1.5,'Fontsize',14,...
    'XTick', []);
xlim([-5,75])
ax.XAxis.Visible='off';
L = legend([a1,a2],{'Preterm at birth','Preterm at TEA'},...
    'Location','North', 'NumColumns',2)
set(L,'Orientation','horizon','Box','off')
ylabel('Coupling index')


% Obtain the number of individuals who crossed the mean line of SFC in
% full-term during development 
[~, idx_l] = find(mPT1<mFT);
l_num = sum(mPT2(idx_l)>mPT1(idx_l))

[~,idx_h] = find(mPT1>mFT);
num = l_num+sum(mPT2(idx_h)<mPT1(idx_h))