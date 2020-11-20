function stop = outfun(x,optimValues,state,varargin)

global NX NY TIMESTEP PARAM PHI LOGPERM_MEAN


stop=false;
persistent PERMITER OBJITER PERMFILE OBJFILE

if  optimValues.iteration < 56

switch state
    case 'iter'
        PERMITER = [PERMITER  x]
        OBJITER = [OBJITER   optimValues.fval]
          % Make updates to plot or guis as needed
                  figure(1)
        switch PARAM

        case 'FULL'
            LOGPERMFULL = x;
            subplot(7,8,  optimValues.iteration + 1)
            imagesc(reshape(LOGPERMFULL,NX,NY) , [5  10]); 
                    
        case 'REDUCED'
            LOGPERMFULL = PHI * x + LOGPERM_MEAN;                          % approximate log-permeability field 
            subplot(7,8, optimValues.iteration + 1)
            imagesc(reshape(LOGPERMFULL,NX,NY) , [5  10]); 
        end 

            figure(2)
            plot( optimValues.iteration, optimValues.fval,'.-','MarkerSize',6); hold on 

    case 'interrupt'

    case 'init'
        PERMITER = x;
        OBJITER = optimValues.fval;
        PERMFILE = ['PERMITER' , num2str(TIMESTEP), '.OUT']
        OBJFILE = ['OBJITER' , num2str(TIMESTEP), '.OUT']
        figure(1)
        switch PARAM

        case 'FULL'
            LOGPERMFULL = x;
            subplot(7,8,1)
            imagesc(reshape(x,NX,NY) , [5  10]); 
                    
        case 'REDUCED'
            LOGPERMFULL = PHI * x + LOGPERM_MEAN;                          % approximate log-permeability field 
            subplot(7,8,1)
            imagesc(reshape(LOGPERMFULL,NX,NY) , [5  10]); 
        end 

            figure(2)
            plot( optimValues.iteration, optimValues.fval,'.-','MarkerSize',6); hold on 

          % Setup for plots or guis
    case 'done'
            FIDPERM = fopen(PERMFILE,'w');
            fprintf(FIDPERM,'%d \n',PERMITER);
            fclose(FIDPERM);

            FIDOBJ = fopen(OBJFILE,'w');
            fprintf(FIDOBJ,'%d \n',OBJITER);
            fclose(FIDOBJ);
            
            close all
          % Cleanup of plots, guis, or final plot
    otherwise
    end

end