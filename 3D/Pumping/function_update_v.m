function out = function_update_v(A,d_obs,alpha);

a = A'*A;
a = a + 0.0001*norm(a)*eye(size(a));
out = pinv(a) * A'*d_obs;

