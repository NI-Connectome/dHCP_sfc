% Brain-behavior analysis
clc; clear; close all

%% Data loading
load ./out/FTPT_cognitive/bsid_data.mat
load ./data/dHCP/FT_436.mat info R
CI_FT = R(:, bsid_FT.idx_FT436);
info_FT = info{bsid_FT.idx_FT436, 3:4}; % [PMA sex]

load ./data/dHCP/PT2_117.mat info R
CI_PT = R(:, bsid_PT2.idx_PT2117);
info_PT = info{bsid_PT2.idx_PT2117, 3:4}; % [PMA sex]

%% Region level
roi = 90;
%% Correlation
for j = 1:3
for i = 1:roi
    %%%%%%%%% FT
    y1 = bsid_FT{:, 2+j};    
    x1 = CI_FT(i,:)';     
    x2 = info_FT;   
    
    % regress out sex and PMA from SFC.
    X = regress_cov(x1,x2);
    
    [rho(i,j), pval(i,j)] = corr(y1, X);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This code block is used to generate data for Figure 4bcde
    if i == 62          % IPL.R
        if j==1
            fig41(:,1) = X;
            fig41(:,2) = y1;
        else if j==2
            fig41(:,3) = X;
            fig41(:,4) = y1;  
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%% PT
    y1 = bsid_PT2{:, 2+j};    
    x1 = CI_PT(i,:)';     
    x2 = info_PT;   
    
    % regress out sex and PMA from CI.
    X = regress_cov(x1,x2);
    [rho(i,j+3), pval(i,j+3)] = corr(y1, X);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This code block is used to generate data for Figure 4bcde
    if i==69 & j==1     % PCL.L
            fig42(:,1) = X;
            fig42(:,2) = y1;
    end
    if i==53 & j==2     % IOG.L
            fig42(:,3) = X;
            fig42(:,4) = y1;  
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
end

% BH correct 
for i = 1:3
    pfdr_a(:,i) = mafdr(pval_a(:,i),'BHFDR',true);
end
