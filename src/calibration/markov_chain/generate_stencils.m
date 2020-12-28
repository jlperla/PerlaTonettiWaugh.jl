function [L_1_minus, L_2] = generate_stencils(z)
    assert(z(end) > 0);
    assert(abs(z(1)) < 1E-10);
    Delta = z(2) - z(1); %Only using a uniform grid for now
    M = length(z);
    M_bar = M+2;

    L_2 = spdiags([ones(M, 1) -2*ones(M, 1) ones(M, 1)], [0 1 2], M, M_bar)/(Delta^2);
    L_1_minus = spdiags([-ones(M, 1) ones(M, 1) zeros(M, 1)], [0 1 2], M, M_bar)/Delta;
end