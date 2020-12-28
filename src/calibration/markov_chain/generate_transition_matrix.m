function [L, f_bar] = generate_transition_matrix(g, mu, upsilon, theta, z, L_1_minus, L_2)
    assert(g - mu > 0);
    assert(upsilon > 0);
    assert(z(end) > 0);
    z_bar = z(end);
    M = length(z);
    M_bar = M+2;
    
    f_bar_pdf = @(z_tilde) theta*exp(-theta*z_tilde)./(1-exp(-theta*z_bar));
    f_bar = f_bar_pdf(z) /sum(f_bar_pdf(z) );
    %%Stationary distribution as the truncated f_bar_pdf
    L_extension = (mu-g) * L_1_minus + upsilon^2/2 * L_2;

    %Start with the interior
    L = L_extension(:, 2:M_bar-1);
    L(M,M-1) = L(M,M-1) + L_extension(M, M_bar); %Manually add in the reflecting barrier at the top
    flux_adopt = L_extension(1,1); %The flow crossing the boundary can be read from the extension
    %L(1,2) = L(1,2) + flux_adopt; %If just reflecting, it would be this
    L(1,:) = L(1,:) + flux_adopt * f_bar; %Instead of reflecting, spreading over the entire distribution in proportion to the empirical pdf.
end