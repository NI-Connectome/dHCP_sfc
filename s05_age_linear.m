
clear all
close all
addpath ./funct

load ./data/NetsInfo.mat
load ./data/dHCP/motions.mat
col = [242, 231, 99;
		218, 52, 38;
		167, 206, 56;
		220, 219, 221;
		34, 129, 128;
		118, 203, 219;
		205, 159, 201;
		58, 66, 143]./255; % color data
legends = {'VIS', 'SM',	'DA',	'VA',	'LM',	'FP',...
    'DMN', 'SUB'};
limy = 0.3; 

%% Dataloading
dataset = 1;
switch dataset
    case 1
        load ./data/dHCP/FT_436.mat
        moti = mo_FT; % motion data
    case 2
        load ./data/dHCP/PT1_120.mat
        moti = mo_PT1;
    case 3
        load ./data/dHCP/PT2_117.mat
        moti = mo_PT2;
end
ga = info.birth_age;
pma = info.scan_age;
sex = info.sex;

%% whole brain
whR = mean(R);
% [r,p] = linear_corr(pma,whR);
y1 = whR';
x1 = pma;
x2 = [sex, moti, ga];       % covariates

y = regress_cov(y1, x2);
[b, bint, r, rint, stats] = regress(y, [ones(size(y)), pma]);


yy = [ones(size(y)), pma]*b;   

correlation_matrix = corrcoef(pma, y); 
corr = correlation_matrix(1, 2); 

loc = [find(pma == min(pma),1), find(pma == max(pma),1)];
x_pma = [pma(loc(1)), pma(loc(2))];
y_pred = [yy(loc(1)), yy(loc(2))];

h = figure(1)
s = scatter(pma,y,'filled');
s.MarkerFaceAlpha = 0.5;
hold on
Line = plot(x_pma,y_pred);

if stats(3)>0.05
    set(Line, 'LineStyle', ':','LineWidth', 2);
else
    set(Line, 'LineWidth', 2);
end

values = ['r=',num2str(corr,'%.3f'),', p=', num2str(stats(3),'%.3f')];
text(max(pma)-2,0.01,values,'fontangle','italic')

fprintf('### linear correlation\nR^2 \tPearson r \tp value \tF value \n');
fprintf('\t\t\t--- whole brain ---\n%.3f \t%.3f \t\t%.3f \t\t%.3f \t%.3f\n', ...
    stats(1), corr, stats(3), stats(2),b(2))

xlabel('PMA(weeks)')
ylabel('Coupling Index')

box off
ylim([0, 0.4])
xlim([min(pma)-0.5, max(pma)+0.5])
set(gca,'TickDir', 'out', 'TickLength', [.01 .01],'FontWeight', 'bold',...
    'fontname','arial','LineWidth',2.0);

%% Net level
clear values
for i=1:8
    idx = find(net==i);
    ntR(i,:) = mean(R(idx,:));
end

fprintf('\t\t\t---  net level  ---\n')
figure
hold on

for i = 1:8
    y1 = ntR(i,:)';
    y = regress_cov(y1, x2);
    
    [b, bint, r, rint, stats] = regress(y, [ones(size(y)), pma]);

    yy = [ones(size(y)), pma]*b;   
    y_pred = [yy(loc(1)), yy(loc(2))];
    Line = plot(x_pma,y_pred);
    

    correlation_matrix = corrcoef(pma,y); 
    netcorr(i,1) = correlation_matrix(1, 2); 

    fprintf('%.3f \t%.3f \t\t%.3f \t\t%.3f \t%.3f \n', ...
        stats(1), corr, stats(3), stats(2),b(2))
    
    netp(i,1) = stats(3);
    
    if stats(3)>0.05
        set(Line, 'LineStyle', ':','LineWidth', 2,  'color', col(i,:));
    else
        set(Line, 'LineWidth', 2, 'color', col(i,:));
    end
     
    values{i,:} = {['r=',num2str(corr,'%.3f')],['p=', num2str(stats(3),'%.3f')]};
    max_y(i,1) = max(y_pred);
    
end

% scatter plot
for i = 8:-1:1
    y1 = ntR(i,:)';
    y = regress_cov(y1, x2);
    s = scatter(pma,y,[],col(i,:),'filled');
    s.MarkerFaceAlpha = 0.2;
    uistack(s,'bottom')
end

[~,idx] = sort(max_y,'descend');


xlabel('PMA(weeks)')
ylabel('Coupling Index')

box off
ylim([0, 0.4])
xlim([min(pma)-0.2, max(pma)+0.2])
set(gca,'TickDir', 'out', 'TickLength', [.01 .01],'FontWeight', 'bold',...
    'fontname','arial','LineWidth',2.0,'Fontsize',16);

%% region level
clear corr
for i = 1:90
    y1 = R(i,:)';
    
    y = regress_cov(y1, x2);
    [b, bint, r, rint, stats] = regress(y, [ones(size(y)), pma]);
    
    
    roip(i,1) = stats(3);

    correlation_matrix = corrcoef(pma, y); 
    roicorr(i,1) = correlation_matrix(1, 2); 

end

%% BH correct
fixedp_net = mafdr(netp,'BHFDR',true);
fixedp_roi = mafdr(roip,'BHFDR',true);