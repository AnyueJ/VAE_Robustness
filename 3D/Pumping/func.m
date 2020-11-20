function [f,g] = func(x)
load logmean logmean
load Dict Dict
KXPAR = x;
kxfull = Dict * x + logmean;
load IND2 IND2
IND = IND2;
% figure 
% colormap jet
% imagesc(reshape(kxfull,[100,100]));
[P_update Grad] = forward_main(kxfull,Dict,KXPAR,IND,1);


load OBS OBS
vec1 = 1/norm(OBS) * ones(length(OBS),1);
Cn = diag(vec1);


misfit = P_update(IND)-OBS;

f = 0.5*(misfit)'*Cn*(misfit);
g = Grad'*Cn*misfit;
end

