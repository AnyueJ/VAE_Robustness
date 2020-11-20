% one phase Pumping Test
%==========================================================================
clear all; close all; clc; fclose('all');
colormap(jet);
MainDir = cd;
startup

%==========================================================================
%%%%%%%%%%%%  RESERVOIR DOMAIN & MODEL SETUP
nx = 100;                                                                   % number of blocks in each direction
ny = 100;
nz = 1; 

load train train
LOGPERM_MEAN = mean(train,2);
[U,S,V] = svd(train-LOGPERM_MEAN);

save LOGPERM_MEAN LOGPERM_MEAN
Dict = [];
Dict = U(:,1:16);
save Dict Dict
load test30 test30
TRUE = test30(:,15);
TRUE(TRUE<5) = -2.3948;
TRUE(TRUE>5) = 0.6009;
figure (1)
imagesc(reshape(TRUE,[100,100]))

%%%%%%%%%%%%  Generating Observations

[P_true] = pressure_calculation(TRUE);

load IND1 IND1
IND = IND1;
OBS_P = P_true(IND);

A_hard = zeros(length(IND),length(TRUE(:)));
for i=1:length(IND)
    A_hard(i,IND(i)) = 1;
end
OBS_HD = A_hard*TRUE(:);

OBS = [OBS_P;OBS_HD];

%%=========================================================================
%%%%%%%%%%  generating Learning data
s_min = min(TRUE);
s_max = max(TRUE);



%%=========================================================================
%  Inversion
ITER_NUM = 0;
k = 16
v=zeros(k,1);
vlast = ones(k,1);
changeratio = norm(Dict*v-Dict*vlast)/norm(Dict*vlast);
Dict0 = Dict;
while(changeratio > 0.01)
    vlast = v;
    ITER_NUM = ITER_NUM + 1
    [d_obs A] = func_linearization(Dict, v, A_hard, OBS,IND);
    v = function_update_v(A,d_obs);
    figure (2)
    colormap jet
    subplot(4,5,ITER_NUM)
    imagesc(reshape(Dict*v+LOGPERM_MEAN,nx,ny)');
misf = sqrt(sum(sum((d_obs-A*v).^2)));
    title({['Ob value = ',num2str(misf)],['True misfit = ',num2str(norm(Dict*v+LOGPERM_MEAN-TRUE))]});
    changeratio = norm(Dict*v-Dict*vlast)/norm(Dict*vlast);
end
PCArecon = Dict*v+LOGPERM_MEAN;
save PCArecon PCArecon


[P_true] = pressure_calculation(PCArecon);


load IND1 IND1
IND = IND1;
OBS_P = P_true(IND);

A_hard = zeros(length(IND),length(TRUE(:)));
for i=1:length(IND)
    A_hard(i,IND(i)) = 1;
end
OBS_HD = A_hard*TRUE(:);

OBS = [OBS_P;OBS_HD];
obsPCA = OBS;
save obsPCA obsPCA