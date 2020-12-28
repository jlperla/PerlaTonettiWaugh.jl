function [zzz, welfare] = compute_growth_fun_cal(params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given parameter values, computes eq. and moments associated with it...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the settings that for finding the solution of the eq. note that
% multiple initial guesses are attempted. Most important is the functions
% command which selects ``selection_gbm_constant_Theta_functions`` which
% are those equations which match up with what is in the paper appendix.

stationary_s.g_guess = [0.025 0.005];
stationary_s.z_hat_guess = [2]; 
stationary_s.Omega_guess = [1];
    
stationary_s.equilibrium_functions = @selection_gbm_constant_Theta_functions; %Swaps out the set of equations used in equilbrium.
stationary_s.debug_output = true; %Displays information on each iteration if debug_output = true.
base_options = optimset('fsolve');
stationary_s.fsolve_options = optimset(base_options, 'TolFun', 1e-12, 'TolX', 1e-12, 'Display', 'off'); %Used when finding the root with fsolve    
stationary_s.stop_on_first_solution = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computes the equillibrium growth rate...

stationary_results = calculate_equilibrium(params, stationary_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next three steps compute the firm moments....
% The static firm moments...

[mass_ex, size_ex] = static_firm_moments(stationary_results, params);

% This function here is another that would need to be replicated for julia.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) compute the transition matrix. These functions are stored in the
% /markov_chain folder

[L, f_bar] = generate_transition_matrix(params.gtarget, params.mu, params.upsilon, params.theta, params.zgridz, params.zgridL1, params.zgridL2);

T = 5; %How far ahead to look
P = expm(L * T); %i.e. solving the ODE for the transition matrix of the finer grid

% Coarsen the grid to four quintiles to match up with data
[P_quintile, ~] = coarsen_to_quintiles_four(P, params.zgridz, f_bar);

% This is the entire transition matrix @ calibrated values, so first row is bottom qunitle
% second row is 2, 
%    0.2051    0.2291    0.2858    0.2799
%    0.2346    0.2464    0.2727    0.2463
%    0.2587    0.2823    0.2883    0.1706
%    0.0861    0.1418    0.3148    0.4573
%
% Then the calibration scheme only targets the bottom two rows!
% If to port to julia, the first step would be to take calibrated values
% and then get julia versions of generate_transition_matrix,
% corsent_to_quintiles to lead to output like above. 

% only grab the top part...
nquants = 2;
upper_quant_vector = P_quintile(end-(nquants-1):end,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output to record stuff

welfare = stationary_results.U_bar;

home_share = stationary_results.lambda_ii;

growth_rate = stationary_results.g;

zzz = [growth_rate,home_share,mass_ex,size_ex,upper_quant_vector(:)']';
