function [obj A] = func_linearization(x, A_hard,OBS,IND);


save x x
system(['python VAErecon.py'])
load VAErecon VAErecon
VAErecon = double(VAErecon);
% figure 
% colormap jet
% imagesc(reshape(kxfull,[100,100]));
[P_update Grad] = forward_main(VAErecon,x,IND,1);  % forward run for jacobian calculation



system(['python VAEjacob.py'])
load VAEjacob VAEjacob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = [Grad;A_hard*VAEjacob];

dd_1 = [P_update(IND);  A_hard * VAErecon'];
obsrecon = dd_1;
save obsrecon obsrecon
d_o = OBS - dd_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   adding Cn

vec1 = 1/norm(OBS(1:end/2)) * ones(length(OBS)/2,1);
vec2 = 1/norm(OBS(end/2+1:end)) * ones(length(OBS)/2,1);
vec = [vec1;vec2];
Cn = diag(vec);

obj = norm(Cn*d_o)

A = -2*(Cn * G)' *d_o;