%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Master Code: Innovation, Firm Dynamics and International Trade, Atkeson and Burstein, February 2010%%%%%%%%%
%% See Online appendix for notes on codes and other details of quantitative model%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1- Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
%MATRIX OF PARAMETER VALUES ; different parameterizations considered in the paper
%b      annualr     nx      shNT       lambda   slopelarge    nf

MATpars=[...


30       0.0283      1.0448     0.4           1      -0.25        0.1


];

%%% Note that I used calibrate_nx to figure out the exporter fixed costs
%%% so trade and firm moments are fit. annualr is from PTW

allcases=[1]; %choose rows to run ... STILL NEED TO CHOOSE newDratio below

%allcases=[3 4 18 19]; %choose rows to run ... STILL NEED TO CHOOSE newDratio below


for indexcase=1:length(allcases)
    
disp(' ')
disp('Computing Atkeson Burstein (2010) model')
    
tic;    

choosecase=allcases(indexcase);   

per=6;                         % Periods per year, equal to 1/DELTA in attached notes
bold=MATpars(choosecase,1);    % Elasticity of process innovation cost (controls 1/elasticity of growth rate -- see attached notes)
annualr=MATpars(choosecase,2); % Annual interest rate; 
lambda=MATpars(choosecase,5);  % Share of labor in production of research good
nx=MATpars(choosecase,3);      % Export fixed cost
anndelta=0.02;               % Annual depreciation rate;
sigma=0.25;                    % Std.Dev of shocks to productivity;
rho=3.17;                         % Elasticity of substitution; PTW

slope=MATpars(choosecase,6);   % Calibrated slope of employment-based distribution for large firms -- choose to match slope for 1000-5000 firms 
shtrade=0.1063;                 % Calibrated share of trade in output; PTW
shNT=MATpars(choosecase,4);    % Calibrated share of employment of exporters; 

%newDratio=0.9995;               % Ratio of new to old D , 0.9995 is the value we use when we consider a small change
newDratio=0.9331;                % "Large change" (maximum to guarantee convergence under b=10)
%newDratio=0.7;                % "Larger change" (maximum to guarantee convergence under b=30)

dotransition=1;                % 1 to compute the transition dynamics  -- note that computing this can take around 4 hours

z0=0;                          % z of entrants; 
L=1;                           % Labor force
n=1;                           % Entry cost 
nf=MATpars(choosecase,7);      % Fixed cost of operation;

T=floor(1000*per^0.5);         % Grid size: 2*T+1 
tol=1e-7;                      % Tolerance in convergence
periodmass=round(5000*per);    % Number of iterations when computing Statinoary Distributions
pertran=round(3000*per^0.5);   % Number of periods in transition dynamics

% Other parameters
r=(1+annualr)^(1/per)-1;       % Interest Rate 
delta=1-(1-anndelta)^(1/per);  % Period probability of exogenous death
beta=1/(1+r);                  % Discount factor

s=sigma/per^0.5;                     % Compute discrete step of progress/regress relative to trend. See attached note.
yy=exp((slope-1)*s);                 % See attached notes on this         
qcal=1/(1-yy^2)*(yy/(1-delta)-yy^2); % Calibrate q for large firms to match slope of firm size distribution

b=bold*2*s*per;                      % Make this adjustment so that elasticity of growth rate to incentives to innovate does not depend on per (see note)

D_1_rho=shtrade/shNT/(1-shtrade/shNT);  % D^(1-rho) chosen to match export share given employment share of exporters

grid=[-(T-1):1:(T-1)]';

vecemployD=exp(z0+s*grid);          % Normalized employment

%disp(' ')
%disp(['Parameter values: primitive b = ',num2str(bold),' , r = ',num2str(annualr),' , lambda = ',num2str(lambda),' , per = ',num2str(per),' , target employ share of exp = ',num2str(shNT),' , target slope large firms = ',num2str(slope)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2- Computing the Initial and final steady States %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('Initial Steady State')

calib=1;
Pars=[per annualr n delta nx z0 T b bold nf rho lambda L beta qcal s];         
H0=0; V0=0; Pid0=0.5;    % Initial guesses for H0, V0 and Pid0 (V0 will not be used in the initial steady state)
[H,ZSS0,ZZSS0,MSS0,YSS0,wSS0,LPSS0,CSS0,XSS0,NTnoDSS0,NNSS0,VxSS0,VSS0,veczSS0,vecqSS0,shareYexpSS0,shareNexpSS0,shareYexphybSS0,PidSS0,FSS0,employ500,StorevarsSS0]=Steady(Pars,grid,periodmass,D_1_rho,calib,H0,V0,Pid0,vecemployD,0,indexcase);  % Initial Steady State

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('New Steady State')

D_1_rho=newDratio^(1-rho)*D_1_rho;
calib=0;
H0=H; V0=VSS0; Pid0=PidSS0;   % Initial guesses for H0, V0 and Pid0 
[H,ZSS1,ZZSS1,MSS1,YSS1,wSS1,LPSS1,CSS1,XSS1,NTnoDSS1,NNSS1,VxSS1,VSS1,veczSS1,vecqSS1,shareYexpSS1,shareNexpSS1,shareYexphybSS1,PidSS1,FSS1,employ500,StorevarsSS1]=Steady(Pars,grid,periodmass,D_1_rho,calib,H0,V0,Pid0,vecemployD,employ500,indexcase);  % New Steady State

elPid=-log(PidSS1/PidSS0)/log(newDratio);
elOut=-log(YSS1/YSS0)/log(newDratio);
elCons=-log(CSS1/CSS0)/log(newDratio);
elLP=-log(LPSS1/LPSS0)/log(newDratio);
elAggProd=-log(ZZSS1/ZZSS0)/(rho-1)/log(newDratio);
elEntry=-log(MSS1/MSS0)/(rho-1)/log(newDratio);
elAverProdwithD=-log(ZSS1/ZSS0)/(rho-1)/log(newDratio);
elAverProdnoD=-(shareNexpSS0*log(NTnoDSS1/NTnoDSS0)+(1-shareNexpSS0)*log(NNSS1/NNSS0))/(rho-1)/log(newDratio);
elExportshare=-log(shareYexpSS1/shareYexpSS0)/log(newDratio);
ratioIndirect=-log(ZZSS1/ZZSS0)/(rho-1)/log(newDratio)/shareYexpSS0-1;

%disp(' ')
%disp(['Steady-state change elasticity: Y = ',num2str(elOut),' , LP = ',num2str(elLP),' , Agg Prod = ',num2str(elAggProd),' , Aver Prod with D exact = ',num2str(elAverProdwithD),' , Aver Prod w/o D = ',num2str(elAverProdnoD),' , Prod Innov = ',num2str(elEntry)])

%save output

groupvars=[bold;annualr;lambda;per;rho;nx;nf;b;H;newDratio;StorevarsSS0;StorevarsSS1;elPid;elOut;elCons;elLP;elAggProd;elEntry;elAverProdwithD;elAverProdnoD;elExportshare;ratioIndirect;toc];
MATstore(1:length(groupvars),indexcase)=groupvars;
dlmwrite('Storeresults',MATstore,'\t');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('Computing Transition Path')

if dotransition==1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3- Computing the Transition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vecemployF=vecemployD*D_1_rho;
kappa=1/lambda^lambda/(1-lambda)^(1-lambda);            % Constant in the research good price index

% I - Saving the Steady State

vecqtran(:,1)=vecqSS0;
Mtran(1)=MSS0;
masstran(1)=sum(MSS0*veczSS0);

% II- Making the Initial Guess. Here we guess Y, Z, C and w. 

Ytranit(1,1)=YSS0;Ytranit(1,2+pertran:1+per+pertran)=YSS1; 
ZZtranit(1,1)=ZZSS0;ZZtranit(1,2+pertran:1+per+pertran)=ZZSS1; 
Ctranit(1,1)=CSS0;Ctranit(1,2+pertran:1+per+pertran)=CSS1; 
wtranit(1,1)=wSS0;wtranit(1,2+pertran:1+per+pertran)=wSS1;

Ytranit(1,2:1+pertran)=[YSS0 YSS0+(YSS1-YSS0)*[1:pertran-1]/(pertran-1)];       % Guess for Y
ZZtranit(1,2:1+pertran)=[ZZSS0 ZZSS0+(ZZSS1-ZZSS0)*[1:pertran-1]/(pertran-1)];  % Guess for Z
Ctranit(1,2:1+pertran)=[CSS0 CSS0+(CSS1-CSS0)*[1:pertran-1]/(pertran-1)];       % Guess for C 
wtranit(1,2:1+pertran)=[wSS0 wSS0+(wSS1-wSS0)*[1:pertran-1]/(pertran-1)];       % Guess for w  

LPtran=zeros(1,2+pertran);LPtran(1,1)=LPSS0;LPtran(1,end)=LPSS1;
Rtran(1,1)=1/beta;Rtran(1,2+pertran:1+per+pertran)=1/beta; 

Ytran(1,1)=YSS0; Ytran(1,2+pertran:1+per+pertran)=YSS1;
ZZtran(1,1)=ZZSS0; ZZtran(1,2+pertran:1+per+pertran)=ZZSS1;
Ctran(1,1)=CSS0; Ctran(1,2+pertran:1+per+pertran)=CSS1;
wtran(1,1)=wSS0; wtran(1,pertran+2:1+per+pertran)=wSS1; 

% III- Setting the parameters for the loop: weights and Tolerance

if newDratio<0.995
    
    if bold==10  |  bold==30 
        if lambda<1
            speedconv=0.015;       % Weights
        elseif lambda==1
            speedconv=0.02;
        end
    else
        speedconv=0.05;
    end    
    toliter=tol*100;           % Setting the Tolerance
    
elseif newDratio>=0.995
     if bold==10
        speedconv=0.02;
    else
        speedconv=0.04;
     end         
     toliter=tol*5;
end

%% Starting the loop. The final Good is used as the numeraire  %%%%%%

diffiter=10;
iter=1;

while diffiter > .5*10^-2

    %% IV- Using the Value function in the last period, and the guesses for Y, Z, C and W, iterate backwards to find the policy functions 

    Vtran=VSS1*kappa*wSS1^lambda;                         % Value at time T, in terms of the final good (P normalized to 1)
      
    for j=1:per
        matVtran(j)=Vtran(T);                                % Stores the values of the last "per" periods, so that they can be used in the free entry condition (entry takes "per" periods)
    end

    for t=pertran+1:-1:2
        
        Ventry=matVtran(per);                              % Value at entry
     
        Rtran(1,t)=1/beta*Ctranit(1,t+1)/Ctranit(1,t);      % Recovering the Interest Rate
        Rcomp=cumprod(Rtran(1,t:t+per-1));
        Rcomptran=Rcomp(end);                               % Annual Interest Rate   
    
        Pid=1/rho*Ytranit(1,t)/ZZtranit(1,t);               % Profits   (Guess)
        Prtranit=kappa*wtranit(1,t)^lambda;                 % Price of Research Good    (Guess)
        
        Prtran=1/Rcomptran*Ventry/n/per;                % Price of Research Good    (Free entry)
        wtran(1,t)=(Prtran/kappa)^(1/lambda);               % Updated wage  (Price Index)
           
        Vup=[Vtran(2:end); Vtran(end)*exp(s)];            % Value next period when moving up
        Vdown=[0; Vtran(1:end-1)];                        % Value next period when moving down                                    
        Vdiff=Vup-Vdown;
        findpos=find(Vdiff>0);   
      
        Vx=max(0,Pid*vecemployD*D_1_rho-Prtranit*nx);      % Value of exporting    
            
        vecqtran(:,t)=zeros(2*T-1,1);         
        vecqtran(findpos,t)=log(1/Rtran(1,t)*(1-delta)/b/Prtranit/H./vecemployD(findpos).*(Vdiff(findpos)))/b;     % Policy function p (foc for p)
     
        vecqtran(find(vecqtran(:,t)<0),t)=0;                % Policy function
        vecqtran(find(vecqtran(:,t)>1),t)=1;               

        c_qtran=H*vecemployD.*exp(vecqtran(:,t)*b);    % Innovation costs
        
        if bold==1111;vecqtran(:,t)=qcal*ones(2*T-1,1);c_qtran=zeros(2*T-1,1);end                % Inelastic process innovation
        
        Vc=Pid.*vecemployD-(nf+c_qtran)*Prtranit+Vx+(1-delta)/Rtran(1,t)*(vecqtran(:,t).*Vup+(1-vecqtran(:,t)).*Vdown);      % Value if active
    
        Vtran=max(0,Vc);                               % Value
        zbar_ac(1,t)=grid(min(find(Vc>=0)));           % Cuttoff to be active
        zbar_ex(1,t)=grid(min(find(Vx>0)));            % Cuttoff to become an exporter

        for j=per:-1:2
            matVtran(j)=matVtran(j-1);
        end
        
        matVtran(1)=Vtran(T); 
    end

    %% V- Using the initial distribution of productivities and the policy
    %% function, iterate forward and solve for the aggregate variables

    vecztran=MSS0*veczSS0;          % Initial distribution (aggregate productivity)

    for t=2:1+pertran
        vecmup=[0; (1-delta)*vecqtran(1:end-1,t-1).*vecztran(1:end-1,1)];               % Mass of firms transitioning up (Positive shock)
        vecmdown=[(1-delta)*(1-vecqtran(2:end,t-1)).*vecztran(2:end,1); (1-delta)*vecqtran(end,t-1).*vecztran(end,1)];   % Mass of firms transitioning down (Negative shock)
    
        if t<1+per
            M_1=Mtran(1);                               % Mas of entrants (decided one year in anticipation)
        elseif t>=1+per
            M_1=Mtran(t-per);
        end
    
        vecztran(:,1)=vecmup+vecmdown;                   
        vecztran(T,1)=vecztran(T,1)+M_1;
        vecztran(1:T+zbar_ac(1,t)-1,1)=0;   
        
        ZZtran(1,t)=sum(vecztran(:,1).*vecemployD)+sum(vecztran(T+zbar_ex(1,t):end,1).*vecemployF(T+zbar_ex(1,t):end,1));       % Aggregate Productivity

        LPtran=speedconv*ZZtran(1,t)*Ytranit(1,t)*(rho/(rho-1)*wtranit(1,t))^-rho+(1-speedconv)*Ytranit(1,t)/ZZtranit(1,t)^(1/(rho-1));  % Labor in the Final Goods Sector
        LMtran=L-LPtran;                                                                % Labor in the research sector
        Xtran=(1-lambda)/lambda*wtranit(1,t)*LMtran;                            
   
        masstran(t)=sum(vecztran(:,1));                                                 % Mass of firms (pay fixed operating cost nf)
        massexptran=sum(vecztran(T+zbar_ex(1,t):end,1));                                     % Mass of exporters (pay fixed cost of exporting ne)
        sharexptran(t)=sum(vecztran(T+zbar_ex(1,t):end,1).*vecemployF(T+zbar_ex(1,t):end,1))/ZZtran(1,t);    

        c_qtran=H*vecemployD.*exp(vecqtran(:,t)*b);                                     % Innovation costs
        if bold==1111;c_qtran=zeros(2*T-1,1);end                                        % Inelastic process innovation
        totalinnovtran=sum(vecztran(:,1).*c_qtran);                                     % Total expenditures in innovation
        restran=LMtran^lambda*Xtran^(1-lambda);                                         % Research good production        
        Mtran(t)=(restran-masstran(t)*nf-massexptran*nx-totalinnovtran)/n/per;          % Mass of entrants (Market clearing for research good)    
        Ytran(1,t)=LPtran*ZZtran(1,t)^(1/(rho-1));                                      % Updated Y
        Ctran(1,t)=Ytran(1,t)-Xtran;                                                    % Updated C
    end
    
    Pid=1/rho*Ytran./ZZtran;               % Profits     
    LPtran=Ytran./ZZtran.^(1/(rho-1));
    masstran(2+pertran:1+per+pertran)=sum(veczSS1(:,1))*MSS1;
    wtran(1+pertran+1:1+pertran+per)=wSS1;
    Pidtran=1/rho*Ytran./ZZtran./(kappa*wtran.^lambda);
    
    %% VI- Update your guess, iterate until convergence
    
    [diffY,indmaxY]=max(abs(log(Ytranit(1,:)./Ytran(1,:))));                            % Distance Y
    diffC=max(abs(log(Ctranit(1,:)./Ctran(1,:))));                                      % Distance C
    [diffw,indmaxw]=max(abs(log(wtranit(1,:)./wtran(1,:))));                            % Distance w
    diffiter=max([diffY diffC diffw]);                                                  % Total distance

    Ytranit=speedconv*Ytran+(1-speedconv)*Ytranit;                                    % New guess for Y                 
    Ctranit=speedconv*Ctran+(1-speedconv)*Ctranit;                                    % New guess for C 
    ZZtranit=speedconv*ZZtran+(1-speedconv)*ZZtranit;                                 % New guess for Z 
    wtranit=speedconv*wtran+(1-speedconv)*wtranit;                                    % New guess for w
    
    % Welfare numbers
    sharexptran(1)=shareYexpSS0;sharexptran(2+pertran:1+per+pertran)=shareYexpSS1;
    welf0=1/(1-beta)*log(Ctran(1));
    welf1=sum(beta.^[0:pertran-1].*log(Ctran(2:1+pertran)))+beta^(pertran)/(1-beta)*log(Ctran(end));
    welfdiff=exp((welf1-welf0)*(1-beta));
    shortlong=(sharexptran(2)-sharexptran(1))/(sharexptran(end)-sharexptran(1));  

    MATstore(length(groupvars)+1:length(groupvars)+5,indexcase)=[-log(welfdiff)/log(newDratio);shortlong;iter;diffiter;toc];
    dlmwrite('Storeresults',MATstore,'\t');
    MATsharexp(:,indexcase)=sharexptran;
    MATZZ(:,indexcase)=ZZtran;
    MATY(:,indexcase)=Ytran;
    MATC(:,indexcase)=Ctran;
    MATPi(:,indexcase)=Pidtran;
    MATM(:,indexcase)=Mtran;
    save restran MATsharexp MATZZ MATY MATC MATPi MATM
    
    direct_effect = -shtrade*(100*log(newDratio));
    
%     disp(['iter = ',num2str(iter),' diff: all = ',num2str(diffiter),' diff: Y = ',num2str(diffY),', w = ',num2str(diffw) ...
%           ,' index: Y = ',num2str(indmaxY),', w = ',num2str(indmaxw),' , welfare = ',num2str(-log(welfdiff)/log(newDratio))])
      %disp(['iter = ',num2str(iter),' diff: all = ',num2str(diffiter),' diff: Y = ',num2str(diffY),', w = ',num2str(diffw)...
       %   ' , welfare elasticity = ',num2str(-log(welfdiff)/log(newDratio)), ' , welfare = ', num2str(100*(welfdiff-1)),...
       %   ' , welfare_direct = ', num2str(direct_effect)])  
    
    iter=iter+1;
    %log(newDratio)
end

disp(' ')
disp('Welfare Gains in Atkeson Burstein (2010) model')
disp(round(100*(welfdiff-1),3))

%disp(['transition dynamics: equivalent variation / log(Dratio) = ',num2str(-log(welfdiff)/log(newDratio)),' , short/long run elast = ',num2str(shortlong)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting the transition %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(indexcase*2)
% Tlen=length(ZZtran);
% linex=(1:1+per+pertran)/per;
% 
% subplot(3,2,1)
% liney=sharexptran;
% plot(linex(1:Tlen),liney(1:Tlen),'-r','LineWidth',2.8)
% axis tight
% title('Panel A: Exports / GDP')
% 
% subplot(3,2,2)
% liney=ZZtran'/ZZtran(1);
% plot(linex(1:Tlen),liney(1:Tlen),'-r','LineWidth',2.8)
% axis tight
% title('Panel B: ZZ')
% 
% subplot(3,2,3)
% liney=wtran'/wtran(1);
% plot(linex(1:Tlen),liney(1:Tlen),'-r','LineWidth',2.8)
% axis tight
% title('Panel C: w')
% 
% subplot(3,2,4)
% liney=Ytran/Ytran(1);
% liney1=Ctran/Ctran(1);
% plot(linex(1:Tlen),liney(1:Tlen),'-r','LineWidth',2.8)
% hold on
% plot(linex(1:Tlen),liney1(1:Tlen),'-b','LineWidth',2)
% hold off
% axis tight
% title('Panel D: Y / Y0 (r) , C/C0 (b)')
% 
% subplot(3,2,5)
% liney=Pidtran'/Pidtran(1);
% plot(linex(1:Tlen),liney(1:Tlen),'-r','LineWidth',2.8)
% axis tight
% title('Panel E: A')
% 
% subplot(3,2,6)
% liney=Mtran'/Mtran(1);
% plot(linex(1:length(liney)),liney,'-r','LineWidth',2.8)
% axis tight
% title('Panel F: M')
% 
% suptitle(['Case ',num2str(indexcase),': Transition Dynamics, b = ',num2str(b)])
% orient tall

end

end

plotforpaper=0;

%first run rows 3, 4, 18 and 19 with newDratio=0.9995, 
%then run this:

% if plotforpaper==1;    
% 
%     numper=100;
%     figure(1)
%     subplot(2,1,1)
%     plot([1:numper*per]/per,-log(MATsharexp(1:numper*per,1)/MATsharexp(1,1))/log(newDratio),'b','LineWidth',2)
%     hold on
%     plot([1:numper*per]/per,-log(MATsharexp(1:numper*per,3)/MATsharexp(1,3))/log(newDratio),'r--','LineWidth',2)
%     hold off
%     title('Panel A: Moderate elasticity of cost of process innovation, b=30')
%     xlabel('Years')
%     ylabel('-\Delta log(Exports/Output)/ \Delta log D')
%     axis([0 100 0 16])
%     legend('Small entering firms','Large entering firms')
% 
%     subplot(2,1,2)
%     plot([1:numper*per]/per,-log(MATsharexp(1:numper*per,2)/MATsharexp(1,2))/log(newDratio),'b','LineWidth',2)
%     hold on
%     plot([1:numper*per]/per,-log(MATsharexp(1:numper*per,4)/MATsharexp(1,4))/log(newDratio),'r--','LineWidth',2)
%     hold off
%     title('Panel B: Low elasticity of cost of process innovation, b=10')
%     xlabel('Years')
%     ylabel('-\Delta log(Exports/Output)/ \Delta log D')
%     axis([0 100 0 16])
% 
%     orient('Tall')
% 
%     %suptitle('Figure 1 : Transition Dynamics of Exports / Output from a Decline in Marginal Trade Costs')
% end