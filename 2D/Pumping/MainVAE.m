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


load test30 test30

test30(test30<5) = -2.3948;
test30(test30>5) = 0.6009;
TRUE = test30(:,15);
figure (1)
imagesc(reshape(TRUE,[100,100])')
%==========================================================================
%%%%%%%%%%%%  Generating Observations

[P_true] = pressure_calculation(TRUE);

OBS_NUM =30;                                                       %%%% NUMBER OF OBSERVATIONS
% % % IND = randperm(nx*ny);
% % % IND = IND(1:OBS_NUM);                                               %%% INDEX OF OBSERVATONS   fix it for comparisson of methods
load IND1 IND1
IND = IND1;
OBS_P = P_true(IND);

A_hard = zeros(length(IND),length(TRUE(:)));
for i=1:length(IND)
    A_hard(i,IND(i)) = 1;
end
OBS_HD = A_hard*TRUE(:);

OBS = [OBS_P;OBS_HD];
obstrue = OBS;
save obstrue obstrue
%%=========================================================================
%%%%%%%%%%  generating Learning data
s_min = min(TRUE);
s_max = max(TRUE);



%%=========================================================================
%  Inversion
ITER_NUM = 0;
vec1 = 1/norm(OBS(1:end/2)) * ones(length(OBS)/2,1);
vec2 = 1/norm(OBS(end/2+1:end)) * ones(length(OBS)/2,1)/15;
vec = [vec1;vec2];
Cn = diag(vec);
x = zeros(1,16);
%load x x
save x x
system(['python VAErecon.py'])
load VAErecon VAErecon
n_z = 16;
m = 0;
v = 0;
t = 1;
alpha = 5e-1;
beta1 = 0.9;
beta2=0.999;
eps = 1e-8;
objhist = zeros(301,1);
xhist = zeros(301,n_z);


reconhist = zeros(301,10000);
reconhist(1,:) = VAErecon;

xhist(1,:) = x;
while t<300
        
        
        figure (2)
        subplot(20,20,t)
        colormap jet
        imagesc(reshape(VAErecon,100,100)');
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


