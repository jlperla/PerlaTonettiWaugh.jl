%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finds the Steady State
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H,Z,ZZ,M,Y,W_P,LP,C,X,NTnoD,NN,Vx,V,vecz,vecq,shareYexp,shareNexp,shareYexphyb,Pid,F,employ500,Storevars] = Steady(Pars,grid,periodmass,D_1_rho,calib,H0,V0,Pid0,vecemployD,employ500,indexcase);

% 1: Recovering Parameters
per=Pars(1); annualr=Pars(2); n=Pars(3); delta=Pars(4); nx=Pars(5); z0=Pars(6); T=Pars(7); 
b=Pars(8) ; bold=Pars(9) ; nf=Pars(10); rho=Pars(11); lambda=Pars(12);  L=Pars(13); beta=Pars(14) ; qcal=Pars(15); s=Pars(16); 

% 2: Find fixed point Pid

warning off

tol=1e-10; %tolerance in convergence
options=optimset('TolFun',tol,'TolX',tol,'Display','off');
Pid=fsolve('Fixed',Pid0,options,Pars,grid,calib,D_1_rho,H0,V0,vecemployD);  % Finds the fixed point Pid
[F,vecq,Vc,Vx,Vtop,c_q,H,V]  = Fixed(Pid,Pars,grid,calib,D_1_rho,H0,V0,vecemployD);

% 3: Computing masses and other distributions (Stationary Distributions)

vecm0=(grid==z0);
vecmhyb0=(grid==z0);
vecz=vecm0;
veczhyb=vecm0;
    
for t=2:periodmass
    vecmup=[0; (1-delta)*vecq(1:end-1).*vecm0(1:end-1)];
    vecmdown=[(1-delta)*(1-vecq(2:end)).*vecm0(2:end); (1-delta)*vecq(end).*vecm0(end)];    
    vecm1=vecmup+vecmdown; 
    vecm1=vecm1.*(Vc>=0);
    vecz=vecz+vecm1;
    vecm0=vecm1; 
    
    vecmup=[0; beta*(1-delta)*vecq(1:end-1).*vecmhyb0(1:end-1)]; %Hybrid productivity (will not be used in computing allocations)
    vecmdown=[beta*(1-delta)*(1-vecq(2:end)).*vecmhyb0(2:end); beta*(1-delta)*vecq(end).*vecmhyb0(end)];    
    vecmhyb1=vecmup+vecmdown;                        
    vecmhyb1=vecmhyb1.*(Vc>=0);
    veczhyb=veczhyb+vecmhyb1;  
    vecmhyb0=vecmhyb1; 
end

% 4: Computes Aggregate Variables

% Normalized employment
    
vecemployF=vecemployD*D_1_rho;
vecemploy=vecemployD+vecemployF;

% Productivity Indices
    
Z=sum(vecz.*(vecemployD+(Vx>0).*vecemployF));       % Total productivity normalized by entering firms
NT=sum(vecz.*(Vx>0).*(vecemployD+vecemployF));      % Total productivity of exporters normalized by entering firms
NN=Z-NT;                                            % Total productivity of non-exporters normalized by entering firms
NTnoD=sum(vecz.*(Vx>0).*vecemployD);                % Total productivity of exporters for domestic sales normalized by entering firms

% Mass of exporters and non-exporters normalized by mass of entering firms.
    
shareNexp=NT/Z;                          %employment share of exporters
shareYexp=shareNexp*D_1_rho/(1+D_1_rho); %exports / intermediate goods output

shareNexphyb=sum(veczhyb.*(Vx>0).*(vecemployD+vecemployF))/sum(veczhyb .*(vecemployD+(Vx>0).*vecemployF)); %hybrid employment share of exporters
shareYexphyb=shareNexphyb*D_1_rho/(1+D_1_rho);                                                             %hybrid export share
    
% Other Variables

Upsilon=sum(vecz.*(nf+(Vx>0).*nx))+sum(vecz.*c_q)+per*n;
ratPie=Pid*Z/Upsilon;
LM=lambda/(lambda+ratPie*(rho-1))*L;  % Research labor
LP=L-LM;                              % Production labor
M=((LM/Upsilon*((1-lambda)/lambda*(rho-1)/rho)^(1-lambda))^(rho-1)*Z^(1-lambda))^(1/(lambda+rho-2)); % Entering firms
W_P=(rho-1)/rho*(Z*M)^(1/(rho-1));
X=LM*W_P*(1-lambda)/lambda;           
ZZ=M*Z;                               % Aggregate productivity
Y=LP*ZZ^(1/(rho-1));                  % Output  
C=Y-X;                                % Consumption  
D=D_1_rho^(1/(1-rho));
Vars=M*(sum(vecz.*(1+(Vx>0))));       %mass of consumed products

%Other statistics

vecemployz=(vecemployD.*(Vx<=0)+vecemploy.*(Vx>0))/Z*LP;   %employment by firm
veccumemploy=cumsum(vecz.*vecemployz)/LP;                  %cdf of employment
gap500=max(find(veccumemploy<0.5));                        %median employment firm

if calib==1
    employ500=vecemployz(gap500);                          %employment of median firm
end

vecemploynum=vecemployz/employ500*500;                     %change in units of employment size

gapmin=min(find(vecz>0));
gap1000=max(find(vecemploynum<1000));
gap5000=max(find(vecemploynum<5000));

gap20000=max(find(vecemploynum<20000));
gap100000=max(find(vecemploynum<100000));
gap200000=max(find(vecemploynum<200000));

slope1000_5000=log((1-veccumemploy(gap5000))/(1-veccumemploy(gap1000)))/log(vecemploynum(gap5000)/vecemploynum(gap1000));
slope20000_100000=log((1-veccumemploy(gap100000))/(1-veccumemploy(gap20000)))/log(vecemploynum(gap100000)/vecemploynum(gap20000));

% figure((indexcase-1)*2+1)
% 
% subplot(2,1,1)
% xx=log(vecemploynum);
% yy=log(1-veccumemploy);
% hold on
% scatter(xx(gapmin:gap200000),yy(gapmin:gap200000),'.')
% hold off
% xlabel('log employment')
% ylabel('log(1-F(employment))')
% if calib==0
%     legend('old SS','new SS')
% end
% title('Employment-based size distribution')
% 
% subplot(2,1,2)
% hold on
% scatter(log(vecemployD(gapmin:gap200000)),vecq(gapmin:gap200000),'.')
% hold off
% xlabel('z')
% ylabel('q')
% title('Process Innovation')

% 6: Displays Results
disp(' ')
disp('Trade/GDP')
disp(round(shareYexp,3))
disp(' ')
%disp(['Pid = ',num2str(Pid),'  , Z = ',num2str(Z),'  , M = ',num2str(M),' ,  ZZ = ',num2str(ZZ),' , last z = ',num2str(vecz(end)),' , last m = ',num2str(sum(vecm0)),'  LP = ',num2str(LP),',  Y = ',num2str(Y),',  C = ',num2str(C),...
%      ' ,  W/P = ',num2str(W_P),' , Cons. varieties = ',num2str(Vars)]) 
%disp(['exp/GDP = ',num2str(shareYexp),' , NT/N = ',num2str(shareNexp),' , slope 1000 - 5000 = ',num2str(slope1000_5000),' , slope large firms = ',num2str(slope20000_100000),' , q100000 = ',num2str(vecq(gap200000)),' , hyb exp/GDP = ',num2str(shareYexphyb)])

Storevars=[D;Pid;Z;M;ZZ;LP;Y;C;W_P;Vars;shareYexp;shareNexp;slope1000_5000;slope20000_100000;vecq(gap200000);shareYexphyb;NT;NN;NTnoD];