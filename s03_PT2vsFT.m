% 3.3 Overgrowth of structure-function coupling in preterm infants
clc; clear; close all

%% data loading
load ./data/NetsInfo.mat
load('./data/dHCP/PT2_117.mat','R')
PT2 = R;
load('./data/dHCP/FT_436.mat','R')
FT = R;

%% Unpaired t-tests

% whole brain & net
whR_PT2 = mean(PT2);
whR_FT = mean(FT);


[h_wh,p_wh] = ttest2(whR_PT2',whR_FT');

% Net level
for i = 1:8
    idx = find(net == i);
    ntR_FT(i, :) = mean(FT(idx, :));
    ntR_PT2(i, :) = mean(PT2(idx, :));
end

for n = 1:8
    x = ntR_FT(n,:);
    y = ntR_PT2(n,:);
    [h_net(n,1),p_net(n,1)] = ttest2(x',y');
end
fixedp_net = mafdr(p_net,'BHFDR',true);
fixedh_net = fixedp_net<0.05;

% ROI
for r = 1:90
    x = FT(r,:);
    y = PT2(r,:);
    [h_roi(r,1),p_roi(r,1)] = ttest2(x',y');
end
fixedp_roi = mafdr(p_roi,'BHFDR',true);
fixedh_roi = fixedp_roi<0.05;


% PT2-FT
diff = mean(PT2,2)-mean(FT,2);
diff = diff.*(fixedp_roi<0.05);
diff(find(diff==0)) = nan;
