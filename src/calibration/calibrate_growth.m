function zzz = calibrate_growth(xxx, params, moments, flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This organizes the paramter values based on what is passed into the
% function

if length(xxx) == 7 % this is the main calibration
    params.d = xxx(1);
    params.theta =xxx(2);
    params.kappa = xxx(3);
    params.chi = 1/xxx(4);
    params.mu = xxx(5);
    params.upsilon = xxx(6);
    params.sigma = xxx(7);
elseif length(xxx) == 4 % This is the Sampson calibration
    params.d = xxx(1);
    params.theta =xxx(2);
    params.kappa = xxx(3);
    params.chi = 1/xxx(4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computes eq and generates all teh relavent moments

Theta_hat = compute_growth_fun_cal(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(xxx) == 7
    
    moments.quantile_moments = moments.quantile_moments';
    moments.quantile_moments = moments.quantile_moments(:);

    Theta = [moments.other_moments'; moments.quantile_moments];

    omega = eye(length(Theta));

    for xxx = 1:length(moments.other_moments)

        omega(xxx,xxx) = 100;
    
    end

    if flag ~= 0    
        zzz = (1/length(Theta)).*(log(Theta)-log(Theta_hat))'*omega*(log(Theta)-log(Theta_hat));
    else
        zzz = [Theta,Theta_hat];
    end
    
elseif length(xxx) == 4 % for sampson
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Theta = [moments.other_moments'];
    Theta_hat = Theta_hat(1:4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    omega = eye(length(Theta));

    for xxx = 1:length(moments.other_moments)

        omega(xxx,xxx) = 100;
    
    end

    if flag ~= 0
    
        zzz = (1/length(Theta)).*(log(Theta)-log(Theta_hat))'*omega*(log(Theta)-log(Theta_hat));

    else
        zzz = [Theta,Theta_hat];
    end
end
 
    
    
    