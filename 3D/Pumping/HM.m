function [recon] = HM(OBS,Dict,A_hard,alpha)
nx = 40;                                                                   % number of blocks in each direction
ny = 120;
nz = 40;
siz = size(Dict);
vsize=siz(2);
v=zeros(vsize,1);
load LOGPERM_TRUE LOGPERM_TRUE
load logmean logmean
TRUE = LOGPERM_TRUE;
ITER_NUM = 0;
load wellloc wellloc
IND = wellloc;
vlast = ones(vsize,1);
changeratio = norm(Dict*v-Dict*vlast)/norm(Dict*vlast+logmean);
while(changeratio > 0.03 & ITER_NUM <3)
    vlast = v;
    ITER_NUM = ITER_NUM + 1
    [d_obs A] = func_linearization(Dict, v, A_hard, OBS,IND);
    v = function_update_v(A,d_obs);
    figure (2)
    colormap jet
    subplot(4,5,ITER_NUM)
    trueplot = flip(permute(reshape(Dict*v+logmean,nx,ny,nz),[3,2,1]),3);
imagesc(trueplot(:,:,1))
misf = sqrt(sum(sum((d_obs-A*v).^2)));
    title({['Ob value = ',num2str(misf)],['True misfit = ',num2str(norm(Dict*v+logmean-TRUE))]});
    changeratio = norm(Dict*v-Dict*vlast)/norm(Dict*vlast+logmean);
end
recon = Dict*v+logmean;
end

