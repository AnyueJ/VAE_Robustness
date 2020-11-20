function [P ADJOINT] = forward_main(Hyd_Con,coeffs,IND,flag_adj);
P =  pressure_calculation_main(Hyd_Con);

tic;
if flag_adj == 1
epsi_adj = 0.000001;  
ADJOINT = zeros(480,length(coeffs));
HC_perturb = zeros(192000,length(coeffs));

delta = epsi_adj*eye(length(coeffs));
x = coeffs+delta;
save x x
system(['python VAErecon.py'])
load VAErecon VAErecon
VAErecon = double(VAErecon);
HC_perturb = VAErecon';
parfor (i=1:length(coeffs),4)
    P1 = pressure_calculation_main(HC_perturb(:,i));
    P1 = P1(IND);
    P2 = P(IND); 
    ADJOINT(:,i) = (P1 - P2)/epsi_adj;

    
end
toc   
end

