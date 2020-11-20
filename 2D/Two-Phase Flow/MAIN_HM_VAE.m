clear; close all; clc;

global NX NY TIMESTEP PARAM PHI PHI_INV LOGPERM_MEAN

GAMMA = 0;                                                                 % prior weight in the obj. fun. 
CINV_OBS = 1;                                                              % obs covariance  
CINV_PRIOR = 1;                                                            % prior cov
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

% INJECTION WELLS
INJWELLSPEC = 'RATE' ; BHPMAX = 100000 * ones(NWINJ,NT);  NPORE = 1;       % water injection rate specified
   % oil production rate specified
PRODWELLSPEC = 'BHP'  ;                                                    % production BHPs specified  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                        % pixel-based permeability estimation
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


load test30 test30
forward(test30(:,15))

load PRODOBS_POW OBS 
obstrue = OBS;
save obstrue obstrue
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

ITER_NUM = 0;
vec1 = 1/norm(OBS(1:end/2)) * ones(length(OBS)/2,1);
vec2 = 1/norm(OBS(end/2+1:end)) * ones(length(OBS)/2,1)/15;
vec = [vec1;vec2];
Cn = diag(vec);
%x = randn(1,16);
x = zeros(1,16);
%load x x
save x x
system(['python VAErecon'])
load VAErecon VAErecon
n_z = 16;
m = 0;
v = 0;
t = 1;
alpha = 1e-1;
beta1 = 0.9;
beta2=0.999;
eps = 1e-8;
objhist = zeros(301,1);
xhist = zeros(301,n_z);


reconhist = zeros(301,10000);
reconhist(1,:) = VAErecon;

xhist(1,:) = x;
while t<200
        
        figure (2)
        subplot(20,10,t)
        colormap jet
        imagesc(reshape(VAErecon,100,100)');
        xlast = x;
        [obj g] = Gradfunc(x,NBLOCK,PARAM,PEQUIL,SWEQUIL,OBS,...
                                        INJWELLSPEC,INJCONTROL,PRODWELLSPEC,PRODCONTROL,NT,NWINJ,NWPROD,...
                                        CINV_OBS,CINV_PRIOR,xlast,GAMMA);
        objhist(t) = obj;
%          if t>1
%               if abs(obj-objhist(t-1))/abs(objhist(t-1))<0.001
%                   
%                  break
%          
%              end
%           end
if mod(t,50)<1
    alpha = alpha*0.8;
end
        g = g';
        t = t + 1
        m = beta1*m+(1-beta1)*g;
        v = beta2*v+(1-beta2)*(g.*g);
        m_hat = m/(1-beta1^t);
        v_hat = v/(1-beta2^t);
        x = x-(alpha*m_hat./(sqrt(v_hat)+eps))';
        
        save x x
        system(['python VAErecon'])
        load VAErecon VAErecon
        
        save t t
        
        reconhist(t,:) = VAErecon;
        save reconhist reconhist
        save objhist objhist
        
        xhist(t,:) = x;
        save xhist xhist
end  
obssave();
%% PLOT ITERATION RESULTS

