function [v,obj] = Inversion(Dict, ini,TRUE)
colormap(jet);
MainDir = cd;
startup
nx = 100;                                                                   % number of blocks in each direction
ny = 100;
nz = 1; 
k = 4;
nfea = 7;
for nnn = 1:nx*ny;
    if TRUE(nnn,1) == 100
        TRUE(nnn,1) = -2.3948;
    else
        TRUE(nnn,1) = 0.6009;
    end
end
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
vec1 = 1/norm(OBS(1:end/2)) * ones(length(OBS)/2,1);
vec2 = 1/norm(OBS(end/2+1:end)) * ones(length(OBS)/2,1)/15;
vec = [vec1;vec2];
Cn = diag(vec);
%%=========================================================================
%%%%%%%%%%  generating Learning data
s_min = min(TRUE);
s_max = max(TRUE);



%%=========================================================================
%  Inversion
ITER_NUM = 0;

v=Dict'*ini;
vlast = ones(length(v),1);
changeratio = norm(Dict*v-Dict*vlast)/norm(Dict*vlast);
Dict0 = Dict;
while(changeratio > 0.01 & ITER_NUM <7)
    vlast = v;
    ITER_NUM = ITER_NUM + 1
    [d_obs A] = func_linearization(Dict, v, A_hard, OBS,IND);
    v = function_update_v(A,d_obs);
    figure (2)
    colormap jet
    subplot(4,5,ITER_NUM)
    imagesc(reshape(Dict*v,nx,ny));
% misf = sqrt(sum(sum((d_obs-A*v).^2)));
%     title({['Ob value = ',num2str(misf)],['True misfit = ',num2str(norm(Dict*v-TRUE))]});
    changeratio = norm(Dict*v-Dict*vlast)/norm(Dict*vlast);
end
full = Dict*v;
[Pt] = pressure_calculation(full);

OBS_NUM =30;                                                       %%%% NUMBER OF OBSERVATIONS
% % % IND = randperm(nx*ny);
% % % IND = IND(1:OBS_NUM);                                               %%% INDEX OF OBSERVATONS   fix it for comparisson of methods
load IND2 IND2
IND = IND2;
OBS_Pt = Pt(IND);

A_hard = zeros(length(IND),length(TRUE(:)));
for i=1:length(IND)
    A_hard(i,IND(i)) = 1;
end
OBS_HDt = A_hard*full(:);

OBSt = [OBS_Pt;OBS_HDt];
obj = norm(OBS-OBSt);
% obj = (Cn*OBS-Cn*OBSt)'*(Cn*OBS-Cn*OBSt);
end

