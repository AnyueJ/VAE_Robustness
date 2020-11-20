
clear; close all; clc;

global NX NY TIMESTEP PARAM PHI PHI_INV LOGPERM_MEAN

GAMMA = 0;                                                                 % prior weight in the obj. fun. 
CINV_OBS = 1;                                                              % obs covariance  
CINV_PRIOR = 1;                                                            % prior covariance 


% RESERVOIR DESCRIPTION
NX = 100  ;   NY = 100  ;   NZ = 1;                                          % number of blocks in each direction
DX = 10  ;   DY = 10  ;   DZ = 10;                                         % size of blocks in each direction
OBS_INDEX = [1:NX:NX*NY     NY:NX:NX*NY];   

POROSITY = 0.2;                                                            % porosity

% RESERVOIR EQUILIBRIUM DATA
PEQUIL = 3000;                                                             % equilibrium pressure (psi)
SWEQUIL = 0.1;                                                             % equilibrium saturation (psi)

% SIMULATION TIME AND TIME STEPS
DT = 90;                                                                   % time step size (days)
NT = 12;                                                                   % # of timesteps

% INPUT WELL DATA
NWINJ = 4;   NWPROD = 1;                                                 % number of production wells
WATDENSITY = 62.0; OILDENSITY = 45.0;                                      % water and oil density
NC = 2;  STREAM = [0 1];                                                   % number of phases and inj. stream 

% INJECTION WELLS
INJWELLSPEC = 'RATE' ; BHPMAX = 100000 * ones(NWINJ,NT);  NPORE = 1;       % water injection rate specified
  % oil production rate specified
PRODWELLSPEC = 'BHP'  ;                                                    % production BHPs specified  

PARAM = 'REDUCED';   NPARAM = 45;                                          % PARAMETERIZED permeability estimation 
NPARAMEACH = [36,9];
GROUPS = [2,5];
%% DATA PREPROCESSING
SIMT = NT * DT;                                                            % total simulation time (days)
NBLOCK = NX * NY * NZ;                                                     % total number of blocks
NWELL = NWINJ + NWPROD;                                                    % total number of wells    
LX = NX * DX;   LY = NY * DY;   LZ = NZ * DZ;                              % reservoir dimensions (m)    
RESVOL = LX * LY * LZ;                                                     % reservoir volume (m^3)      
POREVOL = RESVOL * POROSITY ;                                              % one pore volume (m^3)

TSTART = zeros(1,NT);                                                      % starting time of each simulation (day)
TEND = DT * ones(1,NT);                                                    % end time of each simulation (day)

%% LOAD PERMEABILITY FIELD INFORMATION DATA
% WORK WITH "LOG(PERMEABILITY)" THROUGHOUT THE CODE
load LOGPERMALL LOGPERMALL
NENS = size(LOGPERMALL,2);


%% PARAMETRIZATION OF THE PERMEABILITY FIELD AS SELECTED
switch PARAM
    case 'FULL'                                                            % pixel-based full permeability field
        diary('LOGFILE_FULL.PRT')
        disp('PIXEL-BASED PERMEABILITY FIELD USED FOR ESTIMATION')
        NPARAM = NBLOCK;
        OPTPAR = LOGPERMFULL;
        OPTPAR_PRIOR = LOGPERMFULL;
    case 'REDUCED'                                                         % reduced permeability 
        diary('LOGFILE_KLT.PRT')
        disp('SVD PERMEABILITY FIELD USED FOR ESTIMATION')
        load PHI_ALL PHI_ALL
        load LOGPERM_MEAN LOGPERM_MEAN
        LOGPERMFULL = LOGPERM_MEAN;
        PHI = [];
        j = 1;
        for i = GROUPS
            
            PHI = [PHI PHI_ALL(:,(i-1)*100+1:(i-1)*100+NPARAMEACH(j))];
            j = j+ 1;
        end
        save PHI PHI  
        DM = LOGPERMFULL - LOGPERM_MEAN;                                   % mean removed 
        PERMPAR = inv(PHI'*PHI)*PHI' * DM;                                            % truncated KLT coefficients
        LOGPERMFULL = PHI * PERMPAR + LOGPERM_MEAN;                        % approximated log-permeability field 
        OPTPAR = PERMPAR;                                                  % optimization parameters 
        OPTPAR_PRIOR = PERMPAR;                                            % prior parameters 
end

%% GENERATE OBSERVATIONS VECTOR 


load PRODOBS_POW OBS 


switch INJWELLSPEC
    case 'RATE'                                                            % when inj. rates are specified 
    IRATE = POREVOL * NPORE ./ SIMT ./ 0.15899;                            % rate to flush 1-PORE vol. during simulation
    INJRATE = IRATE * ones(NWINJ,NT)/NWINJ;                                % plan 1 is contant inj. for all time
    
    % .. reset all zero-rates to .1 bpd
    INJRATE(find(INJRATE == 0)) = 1e-1;
    
    % .. initialize the injection bhp as equal to initial res. conditions
    INJBHP(:,1) = PEQUIL * ones(NWINJ,1);
    INJCONTROL = INJRATE;

    case 'BHP'                                                             % when inj. BHPs are specified    
    BHPRISE = 50;
    INJBHPPLAN1 = PEQUIL + BHPRISE * ones(NWPROD,NT);
    INJCONTROL(1:NWINJ,:) = INJBHP;

end

switch PRODWELLSPEC
    case 'RATE'                                                            % when prod. rates are specified   
    PRATE = POREVOL * NPORE ./ SIMT ./ 0.15899;                             % rate to flush NPORE during simulation
        
    PRODRATE = PRATE * ones(NWPROD,NT)/NWPROD;
    PRODRATE(find(PRODRATE == 0)) = 1e-1;

    PRODBHP(:,1) = PEQUIL * ones(NWPROD,1);
    PRODCONTROL = PRODRATE;
   
     case 'BHP'                                                             % when prod. BHPs are specified 

    % .. define pressure drop from initial condition
    BHPDROP = 50;
    PRODBHP = PEQUIL - BHPDROP * ones(NWPROD,NT);
    PRODCONTROL = PRODBHP;
end

% NT =2;
TIMESTEP = NT;

%% CALL OPTIMIZATION CODE FMINUNC WITH BFGS ALGORITHM
OPTIONS = optimset('Display','iter','GradObj','on','LargeScale','off',...
                   'MaxIter',30,...
                   'Tolx',1e-3,'Diagnostics','on','OutputFcn', @OUTFUN,'HessUpdate','bfgs');
               
[PERMSOL,OBJVAL,EXITFLAG,OUTPUT] = fminunc('OBJFUN_ECLIPSE',OPTPAR,OPTIONS,NBLOCK,PARAM,PEQUIL,SWEQUIL,OBS,...
                                        INJWELLSPEC,INJCONTROL,PRODWELLSPEC,PRODCONTROL,NT,NWINJ,NWPROD,...
                                        CINV_OBS,CINV_PRIOR,OPTPAR_PRIOR,GAMMA);
                                    
%% PLOT ITERATION RESULTS

%plotresult(log(5),log(5000))