function [P ADJOINT] = forward_main(Hyd_Con,coeffs,IND,flag_adj);

P =  pressure_calculation_main(Hyd_Con);
ADJOINT = [];

if flag_adj == 1
epsi_adj = 0.000001;  
    
for i=1:1:length(coeffs)
    
    delta = zeros(1,length(coeffs));
    delta(i) = epsi_adj;
    coeff1 = coeffs + delta;
    x = coeff1;
    save x x
    system(['python VAErecon.py'])
    load VAErecon VAErecon
    VAErecon = double(VAErecon);
    HC_perturb = VAErecon;
    P1 = pressure_calculation_main(HC_perturb);
    P1 = P1(IND);
    P2 = P(IND); 
    ADJOINT = [ADJOINT (P1 - P2)/epsi_adj]; 
end   
end

