% Calculate SFC based on communication models
clear all; clc; close all
load ./data/dHCP/Edis_aal.mat
addpath ./funct
load out\contri.mat
%% load data and Eu distance

dis_E = Edis_FT;
for s = 1:size(dis_E,3)
    dis_E(:,:,s) = dis_E(:,:,s)+dis_E(:,:,s)';
end

load ./data/dHCP/FT_436 SC FC

[~,n_ROI,n_subj] = size(SC);

%% shortest path length
L = -log(SC);
for s = 1:n_subj
    l = L(:,:,s);
    spl(:,:,s) = distance_wei(l);
end
spl(isinf(spl))=0;

%% communicability
for s = 1:n_subj
    A = SC(:,:,s);
    D = diag(degrees_und(A)); 
    SCsymm=D^(-1/2)*A*D^(-1/2);
    cmc(:,:,s) = exp(SCsymm);
end

%% regress
for s = 1:n_subj
    for i = 1:n_ROI
        y = FC(:,i,s);
        x1 = dis_E(:,i);
        x2 = spl(:,i,s);
        x3 = cmc(:,i,s);
        X = [ones(size(y)),x1,x2,x3];
        y(i) = [];
        X(i,:) = [];
        X(:,2:4) = zscore(X(:,2:4));   
        [b,~,~,~,stats] = regress(y,X);
        R(i,s) = stats(1); % SFC
        contri(i,:,s) = b(2:4); % contributions for each regression model
    end
end