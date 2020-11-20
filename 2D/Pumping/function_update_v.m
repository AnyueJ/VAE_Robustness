function v = function_update_v(v,d_obs);

alpha = 1e-1;
beta1 = 0.9;
beta2=0.999;
eps = 1e-8
  m = 0;
    v = 0;      
       
        m = beta1*m+(1-beta1)*g
        v = beta2*v+(1-beta2)*(g*g)
        m_hat = m/(1-beta1^t)
        v_hat = v/(1-beta2^t)
        x = x0-alpha*m_hat/(np.sqrt(v_hat)+eps)
        xin = np.zeros((1,n_z))
        xin[0] = x.T  
        xin = model.generator(xin)
        objhist[t] = actvalue(x)
        reconhist[t] = xin
        xhist[t] = x.T
    

