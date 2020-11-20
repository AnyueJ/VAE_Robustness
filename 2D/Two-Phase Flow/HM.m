function [outputArg1,outputArg2] = HM()
global NX NY TIMESTEP PARAM PHI LOGPERM_MEAN

GAMMA = 0;                                                                 % prior weight in the obj. fun. 
CINV_OBS = 1;                                                              % obs covariance  
CINV_PRIOR = 1;                                                            % prior covariance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
NWINJ = 1;   NWPROD = 4;                                                 % number of production wells
WATDENSITY = 62.0; OILDENSITY = 45.0;                                      % water and oil density
NC = 2;  STREAM = [0 1];                                                   % number of phases and inj. stream 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHOOSE WELL SPECIFICATION TYPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INJECTION WELLS
INJWELLSPEC = 'RATE' ; BHPMAX = 100000 * ones(NWINJ,NT);  NPORE = 1;       % water injection rate specified
% 
% INJWELLSPEC = 'BHP'  ;                                                   % injection BHPs specified  

% PRODUCTION WELLS
% PRODWELLSPEC = 'RATE' ; BHPMIN = 1000 * ones(NWINJ,NT);  NPORE = 1;      % oil production rate specified
PRODWELLSPEC = 'BHP'  ;                                                    % production BHPs specified  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHOOSE PERMEABILITY REPRESENTATION FOR ESTIMATION PURPOSES
% PARAM = 'FULL';                                                          % pixel-based permeability estimation
PARAM = 'REDUCED';                                            % PARAMETERIZED permeability estimation 


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


% LOGPERMFULL = LOGPERMALL(:,randi(3500));
% load startpoint startpoint
% LOGPERMFULL = log(startpoint);
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
%         LOGPERM_MEAN = mean(LOGPERMALL,2);
%         LOGPERM_RESID = LOGPERMALL - repmat(LOGPERM_MEAN,1,NENS);
%         [U S V] = svds(LOGPERM_RESID ./ sqrt((NENS-1)) , NPARAM);
%         PHI = U * sqrt(S); 
%         PHI_INV = diag(sqrt(diag(1./S))) *  U'; 
        %load PHI_ALL PHI_ALL
        load LOGPERM_MEAN LOGPERM_MEAN
        LOGPERMFULL = LOGPERM_MEAN;
         %LOGPERMFULL = LOGPERMALL(:,2454);
         PHI = [];
%         j = 1;
%         for i = GROUPS
%             
%             PHI = [PHI PHI_ALL(:,(i-1)*100+1:(i-1)*100+NPARAMEACH(j))];
%             j = j+ 1;
%         end
        load PHI PHI
        DM = LOGPERMFULL - LOGPERM_MEAN;                                   % mean removed 
        PERMPAR = inv(PHI'*PHI)*PHI' * DM;                                            % truncated KLT coefficients
        LOGPERMFULL = PHI * PERMPAR + LOGPERM_MEAN;                        % approximated log-permeability field 
        OPTPAR = PERMPAR;                                                  % optimization parameters 
        OPTPAR_PRIOR = PERMPAR;                                            % prior parameters 
end

%% GENERATE OBSERVATIONS VECTOR 

% OBS = load('OBS135_PWO');
load PRODOBS_POW OBS 
% obsv = 1./sum(OBS.^2,2);
% CINV_OBS = obsv.*ones(NWINJ+2*NWPROD,12);
%% DESIGN PRODUCTION SCENARIO
% FOR RATE CONTROL:  INPUT RELATIVE VALUE OF WELLS INJ/PROD RATE

% FOR BHP CONTROL:   INPUT EXACT VALUE OF PRESSURE DECREASE (-) OR 
%                    INCREASE (+) RELATIVE TO INITIAL PRESSURE

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
%plotresult(log(50),log(5000))

end

