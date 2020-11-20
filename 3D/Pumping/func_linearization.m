function [obj A] = func_linearization(x, A_hard,OBS,IND);

save x x
system(['python VAErecon.py'])
load VAErecon VAErecon
VAErecon = double(VAErecon);

[P_update Grad] = forward_main(VAErecon,x,IND,1);  % forward run for jacobian calculation


G = Grad;

dd_1 = P_update(IND);
obsrecon = dd_1;
save obsrecon obsrecon
d_o = OBS - dd_1;



vec1 = 1/norm(OBS(1:end/2)) * ones(length(OBS)/2,1);
vec2 = 1/norm(OBS(end/2+1:end)) * ones(length(OBS)/2,1)/15;
%vec = [vec1;vec2];
vec = 1/norm(OBS(1:end)) * ones(length(OBS),1);
Cn = diag(vec);

obj = norm(Cn*d_o)

A = -2*(Cn * G)' *d_o;