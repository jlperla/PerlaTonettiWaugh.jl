function [mass_ex, size_ex] = static_firm_moments(stationary_results, params)
% This function computes the BEJK static firm moments.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) The mass of exporters...

mass_ex = (stationary_results.z_hat)^(-params.theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Compute the relative value of an average exporters domestic shipments to a
% average non-exporter's shipments. 

rng(03281978)
z = (1-rand(100000,1)).^(-1./params.theta);

exporter = z > stationary_results.z_hat;

rel_value_ship = (1./z).^(1-params.sigma); % This is the value of domestic shipments 

size_ex = mean(rel_value_ship(exporter==1))./mean(rel_value_ship(exporter~=1));