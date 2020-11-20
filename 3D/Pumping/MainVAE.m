clear all;
close all; clc;
startup
tic;

nx = 40;                                                                   % number of blocks in each direction
ny = 120;
nz = 40;

load testcase testcase

TRUE = testcase;

figure (1)
trueplot = flip(permute(reshape(TRUE,[nx,ny,nz]),[3,2,1]),3);
imagesc(trueplot(:,:,1))

[P_true] = pressure_calculation(TRUE);

                                                     %%%% NUMBER OF OBSERVATIONS
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

s_min = min(TRUE);
s_max = max(TRUE);



%%=========================================================================
%  Inversion
ITER_NUM = 0;
vec1 = 1/norm(OBS(1:end/2)) * ones(length(OBS)/2,1);
vec2 = 1/norm(OBS(end/2+1:end)) * ones(length(OBS)/2,1);
%vec = [vec1;vec2];
vec = 1/norm(OBS(1:end)) * ones(length(OBS),1);
Cn = diag(vec);
n_z = 32;
x = zeros(1,n_z);
%load x x
save x x
system(['python VAErecon.py'])
load VAErecon VAErecon

m = 0;
v = 0;
t = 1;
alpha = 2e-1;
beta1 = 0.9;
beta2=0.999;
eps = 1e-8;
objhist = zeros(301,1);
xhist = zeros(301,n_z);


reconhist = zeros(301,192000);
reconhist(1,:) = VAErecon;

xhist(1,:) = x;
while t<300
        
        
        figure (2)
        subplot(20,20,t)
        colormap jet
        trueplot = flip(permute(reshape(VAErecon,nx,ny,nz),[3,2,1]),3);
        imagesc(trueplot(:,:,1))
        [obj g] = func_linearization(x, A_hard,OBS,IND);
        objhist(t) = obj;
        t = t + 1
        if mod(t,50)<1
        alpha = alpha*0.9;
        end
        g = g+0.01*x';
        m = beta1*m+(1-beta1)*g;
        v = beta2*v+(1-beta2)*(g.*g);
        m_hat = m/(1-beta1^t);
        v_hat = v/(1-beta2^t);
        x = x-(alpha*m_hat./(sqrt(v_hat)+eps))';
        
        save x x
        system(['python VAErecon.py'])
        load VAErecon VAErecon
        
        save t t
        
        reconhist(t,:) = VAErecon;
        save reconhist reconhist
        save objhist objhist
        
        xhist(t,:) = x;
        save xhist xhist
end







