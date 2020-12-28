% $Rev: 3 $
% $Date: 2017-01-06 08:57:38 -0800 (Fri, 06 Jan 2017) $
% $LastChangedBy: jlperla $

function results = calculate_equilibrium(params, settings)
%% Unpack the model parameters
    rho = params.rho;
    sigma = params.sigma;
    n = params.n;
    theta = params.theta;
    gamma = params.gamma;
    d = params.d;
    kappa = params.kappa;
    zeta = params.zeta;
    eta = params.eta;
    Theta = params.Theta;
    chi = params.chi;
    upsilon = params.upsilon;
    mu = params.mu;
    delta = params.delta;

    %Generate the combinations of the guesses.  This uses meshgrid to generate the permutations and then put them into a vector form.
    [g_guess_vec, z_hat_guess_vec, Omega_guess_vec] = meshgrid(settings.g_guess, settings.z_hat_guess, settings.Omega_guess);
    root_guesses = [g_guess_vec(:) z_hat_guess_vec(:) Omega_guess_vec(:)]; %Now this has a list of the permutations of guesses
    
    load_functions = settings.equilibrium_functions; %Can swap out the equilibrum 
    
    load_functions(); %Executes the given set of functions. Must be done after loading of parameters.
    
    %% Calculate the equilibrium growth rate and z_hat
    %Note that z_hat_guesses and g_guess may be an array of starting values, so this will try the permutations of them.
    g_best = NaN; %The NaN will flag whether a candidate has already been found when running through the permutations of paramters
    z_hat_best = NaN;
    Omega_best = NaN;
    root_guess_success = NaN;
    root_count = 0;

    %Function to check if root is interior,
    is_interior = @(g, z_hat, Omega) (g > 0 && Omega > 0 && z_hat >= 1 && ...
        c(g, z_hat, Omega) > 0 && y(g, z_hat, Omega) > 0 && ...
        lambda_ii(g, z_hat, Omega) >= 0 && lambda_ii(g, z_hat, Omega) <= 1 && ...
        sign(1-gamma) * U_bar(g, z_hat, Omega) >= 0 && ... %utility and welfare sign as gamma > 1
        pi_bar_agg(g, z_hat, Omega) > 0 && pi_min(g, z_hat, Omega) > 0 && z_bar(g, z_hat, Omega) > 0 && ...
        x(g, z_hat, Omega) > 0 && r(g, z_hat, Omega) > g && ...
        L_tilde(g, z_hat, Omega) > 0 && L_tilde(g, z_hat, Omega) < 1 && ...
        nu(g, z_hat, Omega) > 0 && a(g, z_hat, Omega) > 0 && b(g, z_hat, Omega) > 0); 
        %S(g, z_hat, Omega) > 0 && ... #Not sure this is working

    root_guess_size = size(root_guesses);
    
    for i = 1:root_guess_size(1) %For every permutation of the parameters, run the root finder and store the results if it is interior.
        root_guess(1) = root_guesses(i,1); %First is the g guess 
        root_guess(2) = root_guesses(i,2); %Second is the z_hat guess
        root_guess(3) = root_guesses(i,3); %Second is the L guess
        
        try %For robustness, this will catch errors thrown by fsolve.
            if(has_selection == true) %Uses different roots depending on whether kappa = 0 or not.
                [rootval, fval, exitflag] = fsolve(@(val) [vm_root(val(1), val(2), val(3)) z_hat_root(val(1), val(2), val(3)), free_entry_root(val(1), val(2), val(3))], [root_guess(1), root_guess(2), root_guess(3)], settings.fsolve_options);
                g_candidate = rootval(1);
                z_hat_candidate = rootval(2);
                Omega_candidate = rootval(3);                   
            else %Has no selection, so just need a root for Omega
                g_candidate = g_value;
                Omega_candidate = Omega_value;
                z_hat_candidate = 1.0;
                fval = 0;
                rootval = 0;
            end
 
        catch err %Errors include the initial value being undefined, etc.  If so, then assume this is not a valid root, flag the results as undefined, and continue.
            if(settings.debug_output == true)
                display(err.message); %Conditionally display the error result.  This doesn't mean that the equilibrium itself has an error as other initial conditions may work.
            end
            rootval = NaN;
            fval = NaN;
            g_candidate = NaN;
            z_hat_candidate = NaN;
            Omega_candidate = NaN;            
        end

        %Testing the candidate and saving the results if it is interior
        g_threshold = .0001; %The threshold for testing if the g are the same.
        imaginary_threshold = .0001; %If the imaginary component is any larger, the root is rejected.
        is_converged_real = norm(fval) < g_threshold  && (abs(max(imag(rootval))) < imaginary_threshold); %Make sure that this didn't find an imaginary root or simply not converge.
        if(is_converged_real && is_interior(g_candidate, z_hat_candidate, Omega_candidate)) %Check if interior
            if(~isnan(g_best) && (norm(g_best - g_candidate) > g_threshold)) %If g has aready been set and this converged to , there was a root already found
                root_count = root_count + 1; %only increments the root count if this is sufficiently different from the old one.
                 if(settings.debug_output == true)
                    display_s('***ERROR: MULTIPLE INTERIOR ROOTS: g = %0.8f and %0.8f.  Choosing the smallest for now', g_best, g_candidate);
                 end
            else
                root_count = max(1, root_count); %If it was close enough, then set or leave the root count as 1.  The max is to make sure we don't overwrite old multiple roots.
            end
            %Unpack the results.  
             if (g_candidate < g_best) || isnan(g_best) %THIS CHOOSES THE SMALLEST INTERIOR ROOT.  LIKELY MISSING A PARAMETER RESTRICTION ON RHO
                g_best = real(g_candidate); %Taking only the real part since we are getting very small imaginary parts on occasion.
                z_hat_best = real(z_hat_candidate); 
                Omega_best = real(Omega_candidate);
                root_guess_success = root_guesses(i);                
             end
        end

    end
    
    %Check if there were any successes, and calculate additional functions if there were.
    if(root_count == 0)
%         if(settings.debug_output == true)
%             display_s('***ERROR: No equilibrium found');
%         end
        results.g = NaN;
        results.z_hat = NaN;
        results.Omega = NaN;
        results.root_guess_success = root_guess_success;
        results.root_count = root_count;
        results.c = NaN;
        results.y = NaN;
        results.lambda_ii = NaN;
        results.U_bar = NaN;
        results.pi_bar_agg = NaN;
        results.x = NaN;
        results.r = NaN;
        %results.S = NaN;
        results.L_tilde = NaN;
        results.z_bar = NaN;
        results.pi_min = NaN;
        results.nu = NaN;
        results.a = NaN;
        results.b = NaN;
        results.d_g_rho = NaN;        
    else

        %Calculate the other equilbrium functions and package the results
        results.g = g_best;
        results.z_hat = z_hat_best;
        results.Omega = Omega_best;
        results.root_guess_success = root_guess_success;
        results.root_count = root_count;
        results.c = c(g_best,z_hat_best,Omega_best);
        results.y = y(g_best,z_hat_best,Omega_best);
        results.lambda_ii = lambda_ii(g_best,z_hat_best,Omega_best);     
        results.U_bar = U_bar(g_best,z_hat_best,Omega_best);
        results.pi_bar_agg = pi_bar_agg(g_best,z_hat_best,Omega_best);
        results.x = x(g_best,z_hat_best,Omega_best);
        results.r = r(g_best,z_hat_best,Omega_best);
        %results.S = S(g_best,z_hat_best,Omega_best);
        results.L_tilde = L_tilde(g_best,z_hat_best,Omega_best);
        results.z_bar = z_bar(g_best,z_hat_best,Omega_best);
        results.pi_min = pi_min(g_best,z_hat_best,Omega_best);
        results.nu = nu(g_best,z_hat_best,Omega_best);    
        results.a = a(g_best,z_hat_best,Omega_best);    
        results.b = b(g_best,z_hat_best,Omega_best);
        results.d_g_rho = sign((1-chi)* pi_bar_agg(g_best,z_hat_best,Omega_best) / Omega_best / pi_min(g_best,z_hat_best,Omega_best)   - 1); %See D.56;
    end
end