
function  [OBJFUN OBJGRAD] = Gradfunc(x,NBLOCK,PARAM,PEQUIL,SWEQUIL,OBS,INJWELLSPEC,...
                                INJCONTROL,PRODWELLSPEC,PRODCONTROL,NT,NWINJ,NWPROD,CINV_OBS,...
                                CINV_PRIOR,xlast,GAMMA)

global PHI LOGPERM_MEAN

NWELL = NWINJ + NWPROD;
NOBS = NWELL + NWPROD;
P(:,1) = PEQUIL * ones(NBLOCK,1);
SW(:,1) = SWEQUIL * ones(NBLOCK,1);
SGAS = zeros(NBLOCK,1);

switch PARAM
        case 'FULL'
            PERM = exp(OPTPAR);
        case 'REDUCED'
            save x x
            system(['python VAErecon'])
            load VAErecon VAErecon
            PERM = exp(VAErecon)';
end

switch INJWELLSPEC
    case 'RATE'
        INJRATE = INJCONTROL;
    case 'BHP'
        INJBHP = INJCONTROL;
end

switch PRODWELLSPEC
    case 'RATE'
        PRODRATE = PRODCONTROL;
    case 'BHP'
        PRODBHP = PRODCONTROL;
end

% .. run simulation for all timesteps
for TI = 1:NT
    disp(['OBTAINING MODEL PREDICTIONS WITH CURRENT PARAMETER ESTIMATES: TIME STEP ' , ...
        num2str(TI) ' OF ' num2str(NT)]);

%% PRINT OUT PERMEABILITY FIELD INPUT FILE FOR ECLIPSE
    FIDPERM = fopen('PERMFIELD.IN','w');
    fprintf(FIDPERM,'PERMX \n');
    % .. PERM is the updated, back-transformed perm vector
    fprintf(FIDPERM,'%d \n',PERM');
    fprintf(FIDPERM,'/');
    fclose(FIDPERM);
    
%% PRINT OUT WELL COMPLETION DATA FOR THE CURRENT TIME STEP

% OPEN INJECTION WELLS FILE  FOR INJ WELLS
    FIDINJWELLS = fopen('INJWELLS.IN' , 'w');
    switch INJWELLSPEC
        case 'RATE'
% PRINT WELL COMPLETION DATA FOR INJRATE CONTROL CASE                
            fprintf(FIDINJWELLS , 'WCONINJE');
            fprintf(FIDINJWELLS , '\n');
            for INJI = 1 : NWINJ
                fprintf(FIDINJWELLS , ['INJ', num2str(INJI),...
                    '    WAT    OPEN   RATE']);
                fprintf(FIDINJWELLS,'   %d          1*        100000  /  \n'...
                    ,INJRATE(INJI,TI));
            end
            fprintf(FIDINJWELLS, '/');
        case 'BHP'
% PRINT WELL COMPLETION DATA FOR INJBHP CONTROL CASE
        fprintf(FIDINJWELLS , 'WCONINJE');
        fprintf(FIDINJWELLS , '\n');
            for INJI = 1 : NWINJ
                
                fprintf(FIDINJWELLS , ['INJ', num2str(INJI),...
                    '    WAT    OPEN   BHP   2*']);
                fprintf(FIDINJWELLS,'  %d       /  \n' ,INJBHP(INJI,TI));            
                 
                fprintf(FIDINJWELLS , ['INJ', num2str(INJI),...
                    '    WAT    OPEN   RATE']);
            end
            fprintf(FIDINJWELLS, '/');
    end
    fclose(FIDINJWELLS);
% OPEN PRODUCTION WELLS FILE FOR PROD WELLS
    FIDPRODWELLS = fopen('PRODWELLS.IN' , 'w');
  
    switch PRODWELLSPEC
        case 'RATE'
% PRINT WELL COMPLETION DATA  PRODRATE CONTROL CASE
            fprintf(FIDPRODWELLS , 'WCONPROD');
            fprintf(FIDPRODWELLS , '\n');
            for PRODI = 1 : NWPROD
                
                fprintf(FIDPRODWELLS , ['PROD', num2str(PRODI),...
                    '    OPEN      LRAT  3*']);
                fprintf(FIDPRODWELLS , '  %d     1*   500  /   \n' , ...
                    PRODRATE(PRODI,TI));
            end
            fprintf(FIDPRODWELLS , '/');
        case 'BHP'
% PRINT WELL COMPLETION DATA FOR PRODBHP CONTROL CASE
            fprintf(FIDPRODWELLS , 'WCONPROD');
            fprintf(FIDPRODWELLS , '\n');
            for PRODI = 1 : NWPROD
                fprintf(FIDPRODWELLS , ['PROD', num2str(PRODI),...
                    '    OPEN      BHP     5*     ']);
                fprintf(FIDPRODWELLS , '    %d     /   \n' , ...
                    PRODBHP(PRODI,TI));
            end
                fprintf(FIDPRODWELLS , '/');
    end
    fclose(FIDPRODWELLS);
    
%% PRINT OUT STATES TO RESTART FILE FOR NEXT RUN
    FIDSTATE = fopen('RESTART.IN','w');
    fprintf(FIDSTATE,'PRESSURE \n');
    fprintf(FIDSTATE,'   %d  \n', P(:,TI));
    fprintf(FIDSTATE,'/  \n');
    fprintf(FIDSTATE,'SWAT \n');
    fprintf(FIDSTATE,'   %d \n', SW(:,TI));
    fprintf(FIDSTATE,'/  \n');
    fprintf(FIDSTATE,'SGAS \n');
    fprintf(FIDSTATE,'   %d \n', SGAS);
    fprintf(FIDSTATE,'/  \n');
    fclose(FIDSTATE);
%% RUN ECLIPSE IN DOS MODE
     [STAT RESULTS] = dos('$e300 reservoir_adj'); clear RESULTS
     
%% READ ECLIPSE STATE FORECAST 
    disp('READING PRESSURE OUTPUTS');
    [RETURNFLAG  TEXTOUT]= dos(['find /n "PRESSURE" ', 'RESERVOIR_ADJ.F0001']);
    STR0 = strfind(TEXTOUT,'['); 
    STR1 = strfind(TEXTOUT,']');
    HEADLINE_PRES = str2num(TEXTOUT(STR0 + 1 : STR1 - 1)); clear STR0 STR1
%     HEADLINE_PRES = 58692;
    P(:,TI+1) = textread('RESERVOIR_ADJ.F0001','%n',NBLOCK,'headerlines',HEADLINE_PRES);
    disp('...DONE');

    disp('READING SATURATION OUTPUTS');
    [RETURNFLAG  TEXTOUT]= dos(['find /n "SWAT" ', 'RESERVOIR_ADJ.F0001']);
    STR0 = strfind(TEXTOUT,'['); 
    STR1 = strfind(TEXTOUT,']');
    HEADLINE_SWAT = str2num(TEXTOUT(STR0 + 1 : STR1 - 1));  clear STR0 STR1
%     HEADLINE_SWAT = 796998;
    SW(:,TI+1) = textread('RESERVOIR_ADJ.F0001','%n',NBLOCK,'headerlines',HEADLINE_SWAT);
    SO(:,TI+1) = 1-SW(:,TI+1);
    disp('...DONE');
    
%% READ WELL FORECAST DATA FROM OUTPUT FILE
    disp('READING WELL OUTPUTS');
    HEADLINE_WELLS = 5;
    WELLDATA = textread('RESERVOIR_ADJ.A0001','%n',1034,'headerlines',HEADLINE_WELLS);
    disp('...DONE');

    disp('READING ADJOINT OUTPUTS');
    [RETURNFLAG  TEXTOUT]= dos(['find /n "''AJGFN  1''       " ', 'RESERVOIR_ADJ.F0001']);
    STR0 = strfind(TEXTOUT,'['); 
    STR1 = strfind(TEXTOUT,']');
    HEADLINE_ADJ = str2num(TEXTOUT(STR0 + 1 : STR1 - 1)); clear STR0 STR1
%     HEADLINE_ADJ = 797589;
    
    FIDADJ = fopen('RESERVOIR_ADJ.F0001');
    ADGFNTEMP = textscan(FIDADJ,'%n','headerlines',HEADLINE_ADJ);
    ADGFN(:,1) = ADGFNTEMP{1}; clear ADGFNTEMP

    TEXTGAP = 3;
    for OBSI = 2:NWINJ+2*NWPROD
    ADGFNTEMP = textscan(FIDADJ,'%n','headerlines',TEXTGAP);
    ADGFN(:,OBSI) = ADGFNTEMP{1};  clear ADGFNTEMP
    end
    fclose(FIDADJ);
    disp('...DONE');
    
    INJBHP_ADJ(:,:,TI) = ADGFN(:,1:NWINJ);
    WOPR_ADJ(:,:,TI) = ADGFN(:, NWINJ + 1 : NWELL );
    WWPR_ADJ(:,:,TI) = ADGFN(:, NWELL + 1 : NOBS ); 
    clear ADGFN
%% EXTRACT FIELD DATA (SCALARS)

%% ACRONYMS : 

%%    FIRST LETTER:  F = FIELD       ;  W = WELL  ;  
%%    SECOND LETTER: O = OIL         ;  W = WATER ;
%%    THIRD LETTER:  P = PRODUCTION  ;  I = INJECTION ;
%%    FOURTH LETTER: R = RATE        ;  T = TOTAL ;
%%    FIFTH LETTER:  C = CUMULATIVE  ; 
FIELDSUM = 10;
    FOPR(TI)  =  WELLDATA(3); FWPR(TI)  =  WELLDATA(4);                    
    FOPT(TI)  =  WELLDATA(5); FWPT(TI)  =  WELLDATA(6);
    FRPV(TI)  =  WELLDATA(7); FWCT(TI)  =  WELLDATA(8);
    FOIP(TI)  =  WELLDATA(9); FWIP(TI)  =  WELLDATA(10);


%% CUMULATIVE FIELD FIELD DATA    
    FOPT_CUM(TI) = sum(FOPT(:)); 
    FWPT_CUM(TI) = sum(FWPT(:));
    FWCT_CUM(TI) = sum(FWCT(:));
    FOE_CUM(TI) = (FOIP(1)-FOIP(TI))/FOIP(1);
    
%% EXTRACT WELL DATA (VECTORS)   
    WOPR = WELLDATA(FIELDSUM + 1 : FIELDSUM + NWELL);         
    WWPR = WELLDATA(FIELDSUM + 1 * NWELL + 1 : FIELDSUM + 2 * NWELL);       
    WWIR = WELLDATA(FIELDSUM +  2 * NWELL + 1 : FIELDSUM + 3 * NWELL);
    
    WOPT  =  WELLDATA(FIELDSUM +  3 * NWELL + 1 : FIELDSUM + 4 * NWELL);
    WWPT  =  WELLDATA(FIELDSUM +  4 * NWELL + 1 : FIELDSUM + 5 * NWELL);
    WWIT  =  WELLDATA(FIELDSUM +  5 * NWELL + 1 : FIELDSUM + 6 * NWELL);
   
    WWCT  =  WELLDATA(FIELDSUM +  6 * NWELL + 1 : FIELDSUM + 7 * NWELL);
    WBHP  =  WELLDATA(FIELDSUM +  7 * NWELL + 1 : FIELDSUM + 8 * NWELL);

    OILPRODRATE(:,TI) = WOPR(NWINJ + 1 : end);
    WATPRODRATE(:,TI) = WWPR(NWINJ + 1 : end);
    WATINJRATE(:,TI) = WWIR(1:NWINJ);

    OILPRODTOT(:,TI) = WOPT(NWINJ + 1 : end);
    WATPRODTOT(:,TI) = WWPT(NWINJ + 1 : end);
    WATINJTOT(:,TI) = WWIT(1:NWINJ);

    WATERCUT(:,TI) = WWPT(NWINJ + 1 : end);
    INJBHP(:,TI) = WBHP(1:NWINJ);
    PRODBHP(:,TI) = WBHP(NWINJ + 1 : end);
end

%% PREDICTED OBSERVATIONS
    if strcmp(INJWELLSPEC , 'RATE') == 1;
        OBS_STATES(1:NWINJ , :) = INJBHP;
    else
        OBS_STATES(1:NWINJ,:) = WATINJRATE;
    end
    
    if strcmp(PRODWELLSPEC , 'RATE') == 1;        
        OBS_STATES(NWINJ + 1 : NWINJ + NWPROD , :) = PRODBHP;
    else
        OBS_STATES(NWINJ + 1 : NWINJ + NWPROD , :) =  OILPRODRATE;
    end
    
    OBS_STATES( NWINJ+NWPROD+1 : NWINJ + 2*NWPROD , : ) = WATPRODRATE;
    
%% FORM THE PARAMETERS TO BE OPTIMIZED
% .. rows are obs. wells (nwinj,nwprod,nwpro), cols are 1..NT
% .. reshape into 1 row; each row segment corresponds to one timestep


CINV_OBS = ones(9,12);
CINV_OBS(1,:) = 1/norm(OBS(1,:));
CINV_OBS(2:5,:) = 1/norm(OBS(2:5,:));
CINV_OBS(6:9,:) = 1/norm(OBS(6:9,:));

CINV_OBS = reshape((CINV_OBS.^2),NOBS*NT,1);
OBS_STATES = reshape(OBS_STATES(:,1:NT),NOBS*NT,1);
OBS = reshape(OBS(:,1:NT),NOBS*NT,1);
% CINV_OBS = reshape(CINV_OBS(:,1:NT),NOBS*NT,1);
OBS_MISMATCH = (OBS_STATES - OBS)' .* (CINV_OBS)' * (OBS_STATES - OBS);
PRIOR_TERM = (x - xlast) * CINV_PRIOR * (x- xlast)';

OBJFUN = sum(OBS_MISMATCH) + GAMMA * sum(PRIOR_TERM);

JACOBIVEC_OBS = reshape([INJBHP_ADJ  WOPR_ADJ  WWPR_ADJ] , NBLOCK , NOBS * NT );
OBJJACOBIVEC_OBS = 2 * JACOBIVEC_OBS .* CINV_OBS' * (OBS_STATES - OBS) ;
OBJGRAD_OBS = sum(OBJJACOBIVEC_OBS , 2);

OBJGRAD_PRIOR = 2 * GAMMA * CINV_PRIOR * (x - xlast);
system(['python VAEjacob'])
load VAEjacob VAEjacob
switch PARAM
    case 'FULL'                                                            % pixel-based full permeability field
        OBJGRAD = OBJGRAD_OBS .* PERM + OBJGRAD_PRIOR;
        disp('PIXEL-BASED GRADIENTS ARE USED')
    case 'REDUCED'                                                         % reduced permeability
        OBJGRAD = (OBJGRAD_OBS .* PERM)' * VAEjacob + 0.0009*x+ OBJGRAD_PRIOR;
        disp('PARAMETERIZATION GRADIENTS ARE USED')
end