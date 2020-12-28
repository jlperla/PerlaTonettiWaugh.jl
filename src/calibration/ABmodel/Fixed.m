%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solves for the Fixed Point in Pied  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,vecq,Vc,Vx,Vtop,c_q,H,V] = Fixed(Pid,Pars,grid,calib,D_1_rho,H,V0,vecemployD);

% 1: Recovering Parameters  
per=Pars(1); annualr=Pars(2); n=Pars(3); delta=Pars(4); nx=Pars(5); z0=Pars(6); T=Pars(7); 
b=Pars(8) ; bold=Pars(9) ; nf=Pars(10); rho=Pars(11); lambda=Pars(12);  L=Pars(13); beta=Pars(14) ; qcal=Pars(15); s=Pars(16); 

Vtop=(1+D_1_rho)*Pid/(1+beta*(1-delta)*((exp(s)-exp(-s))*(1/b-qcal)-exp(-s)));    % Value function of the largest firms. We use it as initial guess in the initial steady state. See attached note.  
                                         
if calib==1    
    H=beta*(1-delta)*Vtop*(exp(s)-exp(-s))/b/exp(qcal*b);      % Calibrated H, initial steady state only
    if H==0;error('H = 0, b too high');end
    V=Vtop.*vecemployD;                                    % Initial Guess final steady state       
elseif calib==0 
    V=V0;                                                      % Initial Guess New steady state 
end

% Value Function Iteration 
Vx=max(Pid*vecemployD*D_1_rho-nx,0);                   % Additional Value for exporters      

tol=1e-10;
diffV=100;                      
iterV=1;   

while diffV>tol    
          
    Vup=[V(2:end,1); V(end)*exp(s)];                        % Value if you move up
    Vdown=[0; V(1:end-1,1)];                                % Value if you move down
    Vdiff=Vup-Vdown;
    findpos=find(Vdiff>0);

    % Policy Function:  max -exp(q*b)*H*+beta*(1-delta)*(q*V(z+s)+(1-q)*V(z-s))  
    % FOC : b*exp(q*b)*H=beta*(1-delta)*(V(z+s)-V(z-s)) 
    
    vecq=zeros(2*T-1,1);         
    vecq(findpos)=log(beta*(1-delta)./(b.*H.*vecemployD(findpos)).*Vdiff(findpos))/b;     % Policy function for q    
    vecq(find(vecq<0))=0;                            
    vecq(find(vecq>1))=1;    
    
    c_q=H*vecemployD.*exp(vecq*b);                          % Innovation costs   
    
    if bold==1111;                                          % Inelastic process innovation            
        vecq=qcal*ones(2*T-1,1);
        c_q=zeros(2*T-1,1);
    end    
    
    Vc=Pid*vecemployD-nf-c_q+Vx+(1-delta)*beta*(vecq.*Vup+(1-vecq).*Vdown);  % Value if active
    Vnew=max(0,Vc);                                                          % Updated Value
                        
    diffV=max(abs(log(V(find(Vnew>0))./Vnew(find(Vnew>0)))));          % Iterating on value functions.
        
    V=Vnew;
    iterV=iterV+1;        
end  
    
F=log(per*n/(V(T)/(1+annualr)));      % Free entry