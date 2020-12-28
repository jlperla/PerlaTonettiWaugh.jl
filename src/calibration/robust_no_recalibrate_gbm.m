clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
growth_r_moments = load('./data/growth_and_r_moments.csv', '-ascii');
firm_moments = load('./data/firm_moments.csv', '-ascii');
bejk_moments = load('./data/bejk_moments.csv', '-ascii');
trade_moments = load('./data/trade_moments.csv', '-ascii');
entry_moments = load('./data/entry_moment.csv', '-ascii');


load('cal_params')
disp('Calbirated values computed on date')
disp(T)
disp('')
disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_param = entry_moments(1,1);
params.zeta = 1.00;
params.gtarget = growth_r_moments(2);
params.rho = growth_r_moments(1) - params.gtarget;

M = 500;
z_bar = 7.0;

addpath('./eq_functions');
addpath('./markov_chain');

% Generating grids and stationary distribution
z = linspace(0, z_bar, M); %The grid, if useful at all.
%The following are invariant as long as M and z are fixed
[L_1_minus, L_2] = generate_stencils(z);

params.zgridL1 = L_1_minus;
params.zgridL2 = L_2;
params.zgridz = z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.gamma = 1.0001;
params.n = 10; % number of countries
params.eta = 0; % denomination of adaption costs
params.Theta = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% moments...
quantile_moments = [firm_moments(1,:); firm_moments(2,:)];

quantile_moments = repmat(1-sum(quantile_moments,2),1,4)/4 + quantile_moments;
% This just ensures they add up to one since, the moments are averaged
% accros years, they don't exactly sum up

moments.quantile_moments = quantile_moments;
moments.other_moments = [params.gtarget,trade_moments(1,1),bejk_moments(1,1),bejk_moments(2,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chi_scale = [0.90, 1, 1.10];
closest_chi = [];

for yyy = 1:length(chi_scale)

scale = 1;

load('cal_params')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.d = new_cal(1);
params.theta = new_cal(2);
params.kappa = new_cal(3);
params.delta = delta_param;

params.sigma = new_cal(7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.chi = 1/ (new_cal(4)*chi_scale(yyy));

params.mu = new_cal(5).*scale;
params.upsilon = new_cal(6).*scale;

disp('Chi Values')
disp(params.chi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[baseline, b_welfare] = compute_growth_fun_cal(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

low_tau = (new_cal(1)-1).*0.90 + 1;
params.d = low_tau;

[counterfact, c_welfare] = compute_growth_fun_cal(params);
params.d = new_cal(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_gain = exp((params.rho).*(c_welfare - b_welfare)) - 1;

record_values = [];
cal_params = [];

record_values = [record_values; baseline(1,1),counterfact(1,1),baseline(1,1)-counterfact(1,1),lambda_gain, params.upsilon];

params.dT = (new_cal(1)-1).*0.90 + 1;
    
    %{'theta', 'kappa', 'chi', 'mu', 'upsilon', 'zeta', 'delta', 'N', 'gamma', 'eta', 'Theta', 'd_0', 'd_T'};
cal_params = [cal_params; baseline(1,1), params.theta, params.kappa, params.chi, params.mu, params.upsilon, params.zeta, params.delta...
    params.n, 1.0, params.eta, params.Theta, params.d, params.dT, params.rho, params.sigma];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_param = params.mu;
upsilon_param = params.upsilon;

scale_values = linspace(.1,2,200);

for zzz = 1:length(scale_values)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    scale = scale_values(zzz);
    
    params.mu = mu_param.*scale;
    params.upsilon = upsilon_param.*scale;
    params.d = new_cal(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [baseline, b_welfare] = compute_growth_fun_cal(params);
    
    %cal_params = [cal_params; baseline(1,1), params.d, params.theta,...
     %   params.kappa, 1/params.chi, params.mu, params.upsilon, params.sigma, params.delta];
    
    params.dT = (new_cal(1)-1).*0.90 + 1;
    
    %{'theta', 'kappa', 'chi', 'mu', 'upsilon', 'zeta', 'delta', 'N', 'gamma', 'eta', 'Theta', 'd_0', 'd_T'};
    cal_params = [cal_params; baseline(1,1), params.theta, params.kappa, params.chi, params.mu, params.upsilon, params.zeta, params.delta...
    params.n, 1.0, params.eta, params.Theta, params.d, params.dT, params.rho, params.sigma];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.d = params.dT;
    
    [counterfact, c_welfare] = compute_growth_fun_cal(params);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    lambda_gain = exp((params.rho).*(c_welfare - b_welfare)) - 1;
       
    record_values = [record_values; baseline(1,1),counterfact(1,1),baseline(1,1)-counterfact(1,1),lambda_gain, params.upsilon];
    
end

%disp(record_values)

filename =  join(['./output/robust/gbm/norecalibrate_values_gbm_',num2str(chi_scale(yyy)),'.mat']);

chi_value = params.chi;

save(filename ,'record_values', 'chi_value')

filename_two =  join(['./output/robust/gbm/param_values_gbm_',num2str(chi_scale(yyy)),'.mat']);

save(filename_two ,'cal_params')

[~,idx] = min(abs(cal_params(:,1)-params.gtarget));
disp('Parameters for given chi which deliver growth closest to baseline growth')
header = {'g', 'theta', 'kappa', 'chi', 'mu', 'upsilon', 'zeta', 'delta', 'N', 'gamma', 'eta', 'Theta', 'd_0', 'd_T'};
disp(header)
disp([cal_params(idx,:)])

filename =  join(['../../parameters/calibration_chi_',num2str(round(params.chi,2)),'.csv']);
header = {'theta', 'kappa', 'chi', 'mu', 'upsilon', 'zeta', 'delta', 'N', 'gamma', 'eta', 'Theta', 'd_0', 'd_T', 'rho', 'sigma'};
writecell([header; num2cell(cal_params(idx,2:end))],filename)

closest_chi = [closest_chi; cal_params(idx,2:end)];

end

save('./output/robust/gbm/closest_chi_params' ,'closest_chi')
%writematrix(closest_chi,'./output/robust/gbm/closest_chi_params.csv')

%
%writecell([header; num2cell(closest_chi)],'../../parameters/closest_chi_params.csv')

rmpath('./eq_functions');
rmpath('./markov_chain');
