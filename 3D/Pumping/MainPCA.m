clear all;
close all; clc;

startup


nx = 40;                                                                   % number of blocks in each direction
ny = 120;
nz = 40;
load PHI_ALL PHI_ALL
PHI = PHI_ALL(:,1:32);
load test test
TRUE = test(:,514);
testcase = TRUE;
save testcase testcase
load logmean logmean

figure (1)
trueplot = flip(permute(reshape(TRUE,[nx,ny,nz]),[3,2,1]),3);
imagesc(trueplot(:,:,1))


[P_true] = pressure_calculation(TRUE);


load wellloc wellloc
IND = wellloc;
OBS_P = P_true(IND);

A_hard = zeros(length(IND),length(TRUE(:)));
for i=1:length(IND)
    A_hard(i,IND(i)) = 1;
end
OBS_HD = A_hard*TRUE(:);

%OBS = [OBS_P;OBS_HD];
OBS = OBS_P;
save OBS OBS

s_min = min(TRUE);
s_max = max(TRUE);




ITER_NUM = 0;
vec1 = 1/norm(OBS(1:end/2)) * ones(length(OBS)/2,1);
vec2 = 1/norm(OBS(end/2+1:end)) * ones(length(OBS)/2,1);
%vec = [vec1;vec2];
vec = 1/norm(OBS(1:end)) * ones(length(OBS),1);
Cn = diag(vec1);
n_z = 32;
x = zeros(1,n_z);
%load x x

m = 0;
v = 0;
t = 1;
alpha = 1e-1;
beta1 = 0.9;
beta2=0.999;
eps = 1e-8;
objhist = zeros(501,1);
xhist = zeros(501,n_z);
PCArecon = x*(PHI')+logmean';

reconhist = zeros(501,192000);
reconhist(1,:) = PCArecon;

xhist(1,:) = x;
while t<10
        
        
        figure (2)
        subplot(4,4,t)
        colormap jet
        trueplot = flip(permute(reshape(PCArecon,nx,ny,nz),[3,2,1]),3);
        imagesc(trueplot(:,:,1))
        [obj,d_obs ,g] = func_linearization(x, A_hard,OBS,IND,PHI,logmean);
        objhist(t) = obj;
        t = t + 1

        x = function_update_v(g,d_obs);
        x = x';
        PCArecon = x*(PHI')+logmean';
        
        save t t
        
        reconhist(t,:) = PCArecon;
        save reconhist reconhist
        save objhist objhist
        
        xhist(t,:) = x;
        save xhist xhist
end







