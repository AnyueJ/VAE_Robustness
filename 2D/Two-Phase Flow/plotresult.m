function [] = plotresult(rmin,rmax)
load PHI PHI
load LOGPERM_MEAN
load LOGPERM_TRUE
load INJ INJ
load PROD PROD
load PERMITER 
OBJITER = load('OBJITER12.OUT');
PERMITER = load('PERMITER12.OUT');
NITER = length(OBJITER);
PERMITER = reshape(PERMITER,length(PHI(1,:)),NITER);

LOGPERMFULL = PHI * PERMITER + repmat(LOGPERM_MEAN,1,NITER);               % approximated log-permeability field
LOGPERMFULL(LOGPERMFULL>rmax) = 11;
LOGPERMFULL(LOGPERMFULL<rmin) = 0;
% figure(1);  set(gcf,'OuterPosition',[20 20 1200 1000],'PaperPositionMode','auto') ;
% 
for j = 1:NITER
    figure(1)
    subplot(5,7,j)
    colormap gray
    imagesc(reshape(LOGPERMFULL(:,j), 100, 100)',[rmin,rmax]);
    axis off
    title(['iter = ',num2str(j-1)])
end
% 
figure(2)
plot(OBJITER(1:end))
ylabel('OBJF');
xlabel('Iter');
title('Objective Function Iterations')
% 
% print -f1 -djpeg90 PERM_ITER.jpg;
% print -f2 -djpeg90 OBJ_ITER.jpg;

figure (3)
set(gcf,'Position',[100 100 500 500])
colormap gray
imagesc(reshape(LOGPERMFULL(:,end), 100, 100)',[rmin,rmax]);
hold on

        scatter(INJ(:,1),INJ(:,2),50,'x','r')
hold on


        scatter(PROD(:,1),PROD(:,2),50,'filled','r')
%legend('Injector', 'Producer')
hold off
title(['OBJVALUE = ',num2str(OBJITER(end))])
figure (4)
set(gcf,'Position',[100 100 500 500])
colormap gray
imagesc(reshape(LOGPERM_TRUE, 100, 100)',[rmin,rmax]);
hold on

        scatter(INJ(:,1),INJ(:,2),50,'x','r')
        


hold on

        scatter(PROD(:,1),PROD(:,2),50,'filled','r')

%
hold off
load netaprox
confidence =activations(netaprox,reshape(LOGPERMFULL(:,NITER), 100, 100)',9);

for i = 1:8
    conf(i) = confidence(1,1,i);
end
figure (5)
set(gcf,'Position',[100 100 500 500])
bar(conf)
xlabel('Scenario')
ylabel('Probability')
title('CNN Prediction')
load net
confidence =activations(net,reshape(LOGPERMFULL(:,NITER), 100, 100)',9);

for i = 1:8
    conf(i) = confidence(1,1,i);
end
figure (55)
set(gcf,'Position',[100 100 500 500])
bar(conf)
xlabel('Scenario')
ylabel('Probability')
title('CNN Prediction')
figure (6)
set(gcf,'Position',[100 100 500 500])
colormap gray
imagesc(reshape(PHI*inv(PHI'*PHI)*PHI'*(LOGPERM_TRUE-LOGPERM_MEAN)+LOGPERM_MEAN, 100, 100)',[rmin,rmax]);

logpermplot = LOGPERMFULL(:,end);
obsplot = [];

    NX = 100  ;   NY = 100  ;   NZ = 1;                                          % number of blocks in each direction
    DX = 10  ;   DY = 10  ;   DZ = 10;                                         % size of blocks in each direction (m)
    
    POROSITY = 0.2;                                                            % porosity
    
    %% READ THE TRUE PERMEABILITY FIELD
    
    PERM = exp(logpermplot);
    %% RESERVOIR EQUILIBRIUM DATA
    PEQUIL = 3000;                                                             % equilibrium pressure (psi)
    SWEQUIL = 0.1;                                                             % equilibrium saturation (psi)
    
    %% SIMULATION TIME AND TIME STEPS
    DT = 90;                                                                   % time step size (days)
    NT = 12;                                                                   % # of timesteps
    
    %% INPUT WELL DATA
    NWINJ = 4;   NWPROD = 1;                                                 % number of production wells
    
    %% WELL SPECIFICATIONS
    
    %% INJECTION WELLS
    INJWELLSPEC = 'RATE' ; BHPMAX = 100000 * ones(NWINJ,NT);  NPORE = 2;       % water injection rate specified
    %
    % INJWELLSPEC = 'BHP'  ;                                                   % injection BHPs specified
    
    %% PRODUCTION WELLS
    % PRODWELLSPEC = 'RATE' ; BHPMIN = 1000 * ones(NWINJ,NT);  NPORE = 1;      % oil production rate specified
    PRODWELLSPEC = 'BHP'  ;                                                    % production BHPs specified
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DATA PERPROCESSING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SIMT = NT * DT;                                                            % total simulation time (days)
    NBLOCK = NX * NY;                                                          % total number of blocks
    NWELL = NWINJ + NWPROD;                                                    % total number of wells
    LX = NX * DX;   LY = NY * DY;   LZ = NZ * DZ;                              % reservoir dimensions (m)
    RESVOL = LX * LY * LZ;                                                     % reservoir volume (m^3)
    POREVOL = RESVOL * POROSITY ;                                              % one pore volume (m^3)
    
    %% SIMULATION AND WELL DATA
    TSTART = zeros(1,NT);                                                      % starting time of each simulation (day)
    TEND = DT * ones(1,NT);                                                    % end time of each simulation (day)
    
    %% PRODUCTION SCENARIO
    
    % FOR RATE CONTROL:  INPUT RELATIVE VALUE OF WELLS INJ/PROD RATE
    
    % FOR BHP CONTROL:   INPUT EXACT VALUE OF PRESSURE DECREASE (-) OR
    %                    INCREASE (+) RELATIVE TO INITIAL PRESSURE
    
    switch INJWELLSPEC
        case 'RATE'                                                            % when inj. rates are specified
            IRATE = POREVOL * NPORE ./ SIMT ./ 0.15899;                              % rate to flush 1-PORE vol. during simulation
            INJRATE = IRATE * ones(NWINJ,NT)/NWINJ;                                      % plan 1 is contant inj. for all time
            
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
    
    %% RESERVOIR EQUILIBRIUM/INITIAL CONDITION
    
    P(:,1) = PEQUIL * ones(NBLOCK,1);
    SW(:,1) = SWEQUIL * ones(NBLOCK,1);
    SO(:,1) = (1-SWEQUIL) * ones(NBLOCK,1);
    SGAS = zeros(NBLOCK,1);
    
    for TI = 1:NT
        disp(['TIME STEP ' ,  num2str(TI) , ' CORRESPONDING TO ', ...
            num2str(3*TI), ' MONTHS']);
        
        %% PRINT OUT PERMEABILITY FIELD
        FIDPERM = fopen('TRUE_RESERVOIR.IN','w');
        fprintf(FIDPERM,'PERMX \n');
        fprintf(FIDPERM,'%d \n',PERM');
        fprintf(FIDPERM,'/');
        fclose(FIDPERM);
        
        %% PRINT OUT WELL COMPLETION DATA FOR THE CURRENT TIME STEP
        
        %% OPEN INJECTION WELLS FILE
        FIDINJWELLS = fopen('INJWELLS.IN' , 'w');
        switch INJWELLSPEC
            case 'RATE'
                %% PRINT WELL COMPLETION DATA
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
                %% PRINT WELL COMPLETION DATA
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
        
        %% OPEN PRODUCTION WELLS FILE
        FIDPRODWELLS = fopen('PRODWELLS.IN' , 'w');
        
        %% PRINT WELL COMPLETION DATA
        switch PRODWELLSPEC
            case 'RATE'
                %% PRINT WELL COMPLETION DATA
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
                
                %% PRINT WELL COMPLETION DATA
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
        
        %% PRINT OUT STATES TO RESTART FOR NEXT RUN
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
        [STATUS RESULTS] = dos('$e300 TRUE_RESERVOIR');
        
        %% READ ECLIPSE STATE FORECAST
        [RETURNFLAG  TEXTOUT]= dos(['find /n "PRESSURE" ', 'TRUE_RESERVOIR.F0001']);
        STR0 = strfind(TEXTOUT,'[');
        STR1 = strfind(TEXTOUT,']');
        HEADLINE_PRES = str2num(TEXTOUT(STR0 + 1 : STR1 - 1)); clear STR0 STR1
        disp('READING STATE OUTPUTS');
        %     HEADLINE_PRES = 58692;
        P(:,TI+1) = textread('TRUE_RESERVOIR.F0001','%n',NBLOCK,'headerlines',HEADLINE_PRES);
        
        [RETURNFLAG  TEXTOUT]= dos(['find /n "SWAT" ', 'TRUE_RESERVOIR.F0001']);
        STR0 = strfind(TEXTOUT,'[');
        STR1 = strfind(TEXTOUT,']');
        HEADLINE_SWAT = str2num(TEXTOUT(STR0 + 1 : STR1 - 1));  clear STR0 STR1
        
        %     HEADLINE_SWAT = 796998;
        SW(:,TI+1) = textread('TRUE_RESERVOIR.F0001','%n',NBLOCK,'headerlines',HEADLINE_SWAT);
        SO(:,TI+1) = 1-SW(:,TI+1);
        
        disp('...COMPLETED');
        
        %% READ WELL FORECAST DATA FROM OUTPUT FILE
        [RETURNFLAG  TEXTOUT]= dos(['find /n "PARAMS" ', 'TRUE_RESERVOIR.A0001']);
        STR0 = strfind(TEXTOUT,'[');
        STR1 = strfind(TEXTOUT,']');
        HEADLINE_WELLS = str2num(TEXTOUT(STR0 + 1 : STR1 - 1)); clear STR0 STR1
        
        disp('READING WELL OUTPUTS');
        %     HEADLINE_WELLS = 5;
        WELLDATA = textread('TRUE_RESERVOIR.A0001','%n',730,'headerlines',HEADLINE_WELLS);
        disp('...COMPLETED');
        
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
    
    switch INJWELLSPEC
        
        case 'RATE'
            switch PRODWELLSPEC
                case 'RATE'
                    save TRUTH45_RATE_RATE.mat P SO INJBHP WATINJRATE WATINJTOT ...
                        PRODBHP OILPRODRATE WATPRODRATE OILPRODTOT WATPRODTOT ...
                        WATERCUT FOPR FWPR FOPT FWPT FWCT FOIP FWIP
                case 'BHP'
                    save TRUTH45_RATE_BHP.mat P SO INJBHP WATINJRATE WATINJTOT ...
                        PRODBHP OILPRODRATE WATPRODRATE OILPRODTOT WATPRODTOT ...
                        WATERCUT FOPR FWPR FOPT FWPT FWCT FOIP FWIP
            end
            
        case 'BHP'
            switch PRODWELLSPEC
                case 'BHP'
                    save TRUTH45_BHP_BHP.mat P SO INJBHP WATINJRATE WATINJTOT ...
                        PRODBHP OILPRODRATE WATPRODRATE OILPRODTOT WATPRODTOT ...
                        WATERCUT FOPR FWPR FOPT FWPT FWCT FOIP FWIP
                case 'RATE'
                    save TRUTH45_BHP_RATE.mat P SO INJBHP WATINJRATE WATINJTOT ...
                        PRODBHP OILPRODRATE WATPRODRATE OILPRODTOT WATPRODTOT ...
                        WATERCUT FOPR FWPR FOPT FWPT FWCT FOIP FWIP
            end
    end
    
    set(0,'defaultaxesfontsize',8);
    
    figure; set(gcf,'OuterPosition',[20  50 1300 3*275],'PaperPositionMode','auto')
    
    
    subplot(3,5,1)
    imagesc(reshape(PERM,NX,NY)');
    title('Permeability')
    
    subplot(3,5,2)
    plot([0:3:36], [0 OILPRODTOT])
    title('WOPT')
    
    subplot(3,5,3)
    plot([0:3:36], [0 WATPRODTOT])
    title('FWPT')
    
    subplot(3,5,4)
    plot([3:3:36], mean(INJBHP,1))
    title('Mean INJBHP')
    
    
    IPRINT = [1 2 3 7 13];
    
    for j = 1:length(IPRINT)
        subplot(3,5,5+j)
        Pprint = reshape(P(:,IPRINT(j)),NX,NY)'; colorbar
        
        imagesc(Pprint); % colorbar('southoutside')
        %     title([num2str(3*(IPRINT(j)-1))  ' month'])
        
        subplot(3,5,10+j)
        SOprint = reshape(SO(:,IPRINT(j)),NX,NY)'; colorbar
        imagesc(SOprint,[0 1]);
        %     title([num2str(3*(IPRINT(j)-1))  ' month'])
        
    end
    
    print -f1 -djpeg90 RESULTS_TRUE
    obsnew = [INJBHP ; OILPRODRATE ; WATPRODRATE];
    obsplot = [obsplot obsnew];

    load PRODOBS_POW OBS
    allo = [OBS obsplot];
    minvalue1 = min(allo(1:4,:),[],'all');
    minvalue2 = min(allo(5,:),[],'all');
    minvalue3 = min(allo(6,:),[],'all');
    maxvalue1 = max(allo(1:4,:),[],'all');
    maxvalue2 = max(allo(5,:),[],'all');
    maxvalue3 = max(allo(6,:),[],'all');
    for j = 1:4
        figure (10)
        set(gcf,'Position',[100 100 1500 600])
        subplot(2,4,j)
        
        plot(1:12,obsplot(j,1:12),'LineWidth',2)
        axis([1 12 minvalue1 maxvalue1])
        hold on
        scatter(1:12,OBS(j,:),'filled','g')
        axis([1 12 minvalue1 maxvalue1])
        title(['BHP at INJ', num2str(j)]);
        xlabel('Time Step (90 days)')
        ylabel('Pressure (psi)')
        hold off
        %legend('Reconstruction','Observation','location','best')
    end
    
        figure (10)
        subplot(2,4,5)
        
        
        plot(1:12,obsplot(5,1:12),'LineWidth',2)
        hold on
        scatter(1:12,OBS(5,:),'filled','g')
        title(['Oil rate at PROD1']);
        xlabel('Time Step (90 days)')
        ylabel('Rate (stb/day)')
        hold off
        
        subplot(2,4,6)
        
        plot(1:12,obsplot(6,1:12),'LineWidth',2)
        hold on
        scatter(1:12,OBS(6,:),'filled','g')
        title(['Water rate at PROD1']);
        xlabel('Time Step (90 days)')
        ylabel('Rate (stb/day)')
        hold off
end

