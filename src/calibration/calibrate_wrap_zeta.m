%clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(03281978)
addpath('./eq_functions');
addpath('./markov_chain');

growth_r_moments = load('./data/growth_and_r_moments.csv', '-ascii');
firm_moments = load('./data/firm_moments.csv', '-ascii');
bejk_moments = load('./data/bejk_moments.csv', '-ascii');
trade_moments = load('./data/trade_moments.csv', '-ascii');
entry_moments = load('./data/entry_moment.csv', '-ascii');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('')
disp('')
disp('Calbirated values computed on date')
disp(today('datetime'))
disp('')
disp('')
disp('Calibration with different zetas')
disp('')
disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.gamma = 1.0001;
params.n = 10; % number of countries
params.eta = 0; % denomination of adaption costs
params.Theta = 1;
params.delta = entry_moments(1,1);
params.gtarget = growth_r_moments(2);
params.rho = growth_r_moments(1) - params.gtarget;

M = 500;
z_bar = 7.0;

% Generating grids and stationary distribution
z = linspace(0, z_bar, M); %The grid, if useful at all.
%The following are invariant as long as M and z are fixed
[L_1_minus, L_2] = generate_stencils(z);

params.zgridL1 = L_1_minus;
params.zgridL2 = L_2;
params.zgridz = z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% moments...
quantile_moments = [firm_moments(1,:); firm_moments(2,:)];

quantile_moments = repmat(1-sum(quantile_moments,2),1,4)/4 + quantile_moments;
% This just ensures they add up to one since, the moments are averaged
% accros years, they don't exactly sum up

moments.quantile_moments = quantile_moments;

moments.other_moments = [params.gtarget,trade_moments(1,1),bejk_moments(1,1),bejk_moments(2,1)];
% productivity growth, home share, frac exporters, relative size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeta_scale = [0.75, 1.25];

for yyy = 1:length(zeta_scale)

disp('Zeta Value')

disp(zeta_scale(yyy))

params.zeta = zeta_scale(yyy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.gamma = 1.0001;
initial_val = [3.0220    4.9898    0.1042.*params.zeta    7.8833   -0.0311    0.0483    3.1673];

options = optimset('Display','final','MaxFunEvals',5e4,'MaxIter',1e5);

tic
new_cal =fminsearch(@(xxx) calibrate_growth(xxx,params,moments,1),initial_val,options);
toc

disp('Parameter Values')
disp('d, theta, kappa, 1/chi, mu, upsilon, sigma, delta, rho, zeta')
disp([new_cal,params.delta,params.rho, zeta_scale(yyy)])
disp('')
disp('')
disp('zeta/kappa')
disp(params.zeta/new_cal(3))

all_stuff = calibrate_growth(new_cal,params,moments,0);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('')
disp('Moments: Model Alt Zetas')
disp('Calibration Targets...')
disp('')
disp('Real Rate and Productivity Growth')
disp([params.rho + all_stuff(1,2), all_stuff(1,2)])
disp('')
disp('BEJK Exporter Moments: Fraction of Exporters, Relative Size')
disp([all_stuff(3,2), all_stuff(4,2)])
disp('')
disp('Home Trade Share')
disp([all_stuff(2,2)])
disp('')
disp('Entry Moment')
disp(params.delta)
disp('Firm Moments')
disp(reshape(all_stuff(5:end,2)',2,4))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.d = new_cal(1);
params.theta = new_cal(2);
params.kappa = new_cal(3);
params.chi = 1/new_cal(4);
params.mu = new_cal(5);
params.upsilon = new_cal(6);
params.sigma = new_cal(7);
params.gamma = 1.0;
params.dT = (new_cal(1)-1).*0.90 + 1;


header = {'theta', 'kappa', 'chi', 'mu', 'upsilon', 'zeta', 'delta', 'N', 'gamma', 'eta', 'Theta', 'd_0', 'd_T', 'rho', 'sigma'};

final_cal = [params.theta, params.kappa, params.chi, params.mu, params.upsilon, params.zeta, params.delta...
    params.n, params.gamma, params.eta, params.Theta,params.d, params.dT, params.rho, params.sigma];

filename =  join(['../../parameters/calibration_zeta_',num2str(params.zeta),'.csv']);

writecell([header; num2cell(final_cal)],filename)

end

rmpath('./eq_functions');
rmpath('./markov_chain');












