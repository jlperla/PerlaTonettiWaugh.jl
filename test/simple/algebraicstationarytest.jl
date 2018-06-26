# Unit testing for the algebraic stationary solution, simple model. 

# Arnav Sood: Jun. 26, 2018

using Base.Test
include("..\\..\\src\\simple\\algebraicstationary.jl"); 

# Tests for default parameters 
@test checkparams(simpleparams()) == true # Check that the default parameters are valid.  
@test g(simpleparams()) ≈ 0.021118282603 # Check that the default parameter g value doesn't change. 
@test υ(g(simpleparams()), simpleparams()) ≈ 1.753699551 # Check that the default parameter υ value doesn't change. 
@test v(2, υ(g(simpleparams()), simpleparams()), simpleparams()) ≈ 165.3158126 # Check that the default parameter v value (at state 2) doesn't change. 

# Check parameter validation

# Check g

# Check υ

# Check v (function return)

# Check v (value return)

# Check stationary_algebraic (no z)

# Check stationary_algebraic (with z)