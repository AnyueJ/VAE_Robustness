function [W,group_norm] = WEIGHT_Func(KXPAR,NUM,ITER,WW);

p =   2;
len = length(KXPAR);
epsi = 2^-(ITER+2);
%epsi = 10^-2;
LEN = length(NUM);
group_norm = 0;


for i=1:LEN
    if i==1
        tt = 0;
    else 
        tt = sum(NUM(1:i-1));
    end
    group_size = sum(NUM(1:i))-tt;
        
    vec = KXPAR(tt+1 : tt+group_size);
    WEIGHT = diag(WW{i});
    
    group_norm = group_norm + sum(WEIGHT .* abs(vec).^p )^(1/p);
    
    for j=1:group_size
        p = p;
        
        temp = WEIGHT(j) / (( sum(WEIGHT .* abs(vec).^p))^(1-1/p)+ epsi);
        W(tt + j,tt + j) = temp;
        
    end

end
