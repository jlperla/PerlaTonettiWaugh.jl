function [P_quintile, quantile_indices, quantile_values] = coarsen_to_quintiles_four(P, z, f_bar)
    M = length(z);
    tol = 1e-12; % Based on running 'test_transition_matrix', whose diffs bottom out in magnitude around 1e-11. 
    % Sum over f_bar. 
    f_bar_cdf = cumsum(f_bar);
    % Define the analytical "quantile function." 
    get_quantile = @(x) z(find(f_bar_cdf - x>= -tol, 1, 'first')); % Gives us the grid point with highest cumulative mass <= x. 
    get_quantile_loc = @(x) find(f_bar_cdf - x >= -tol, 1, 'first');
    % Get the "exact" quantile boundaries. 
    first_loc = get_quantile_loc(0.25); 
    second_loc = get_quantile_loc(0.50); 
    third_loc = get_quantile_loc(0.75); 
    fourth_loc = get_quantile_loc(1.0);   
    % Package them.
    first_quintile = 1:first_loc; 
    second_quintile = first_loc + 1: second_loc; 
    third_quintile = second_loc + 1: third_loc; 
    fourth_quintile = third_loc + 1: fourth_loc; 
    % Return exact quintile function to check. 
    quantile_indices = [first_loc, second_loc, third_loc];
    quantile_values = arrayfun(@(x) get_quantile(x), [0.25, 0.50, 0.75]);
    %% Check the quintiles are different
    assert(first_quintile(1) <= first_quintile(end) < second_quintile(1) <= second_quintile(end) <   ...
        third_quintile(1) <= third_quintile(end) < fourth_quintile(1) <= fourth_quintile(end), 'degenerate quintiles'); ...
        assert(first_quintile(1) < first_quintile(end) < second_quintile(1) < second_quintile(end) < third_quintile(1) < ...
        third_quintile(end) < fourth_quintile(1) < fourth_quintile(end), 'degenerate quintiles');
    %% Collapse along columns
    sum_P = zeros(M, 4); 
    for i = 1:M
        sum_P(i, 1) = sum(P(i, first_quintile));
        sum_P(i, 2) = sum(P(i, second_quintile));
        sum_P(i, 3) = sum(P(i, third_quintile));
        sum_P(i, 4) = sum(P(i, fourth_quintile)); 
    end 
    %% Collapse along rows. 
    P_quintile = zeros(4, 4); 
    weights = zeros(1, M);
    % Define the weights. 
    for i = first_quintile
        weights(i) = f_bar(i)/sum(f_bar(first_quintile));
    end 
    for i = second_quintile
        weights(i) = f_bar(i)/sum(f_bar(second_quintile)); 
    end
    for i = third_quintile
        weights(i) = f_bar(i)/sum(f_bar(third_quintile)); 
    end 
    for i = fourth_quintile
        weights(i) = f_bar(i)/sum(f_bar(fourth_quintile)); 
    end 
    
    % Collapse over columns, weighting by the conditional probability
    for j = 1:4
        P_quintile(1, j) = dot(weights(first_quintile), sum_P(first_quintile, j)); 
    end
    for j = 1:4
        P_quintile(2, j) = dot(weights(second_quintile), sum_P(second_quintile, j)); 
    end
    for j = 1:4
        P_quintile(3, j) = dot(weights(third_quintile), sum_P(third_quintile, j)); 
    end
    for j = 1:4
        P_quintile(4, j) = dot(weights(fourth_quintile), sum_P(fourth_quintile, j)); 
    end

end