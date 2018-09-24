# setup, unit circles
lb = 0.0
ub = 2.0
compactifier = Compactifier(lb, ub)
decompactifier = Decompactifier(compactifier)

# values for try
value_before_compactified_1 = lb - 12.1     # out of bounds
value_before_compactified_2 = ub + 19.2     # out of bounds 
value_before_compactified_3 = lb            # lb
value_before_compactified_4 = ub            # ub
value_before_compactified_5 = (lb + ub) / 2 # within bounds
value_before_compactified_6 = lb - 524.1     # out of bounds, more extreme values
value_before_compactified_7 = ub + 595.2     # out of bounds, more extreme values

@test compactifier(-999999999.0) ≈ lb
@test compactifier(999999999.0) ≈ ub
@test decompactifier(lb) ≈ -Inf
@test decompactifier(ub) ≈ Inf
@test compactifier(0.0) ≈ (lb + ub) / 2

@test decompactifier(compactifier(value_before_compactified_1)) ≈ value_before_compactified_1
@test decompactifier(compactifier(value_before_compactified_2)) ≈ value_before_compactified_2
@test decompactifier(compactifier(value_before_compactified_3)) ≈ value_before_compactified_3
@test decompactifier(compactifier(value_before_compactified_4)) ≈ value_before_compactified_4
@test decompactifier(compactifier(value_before_compactified_5)) ≈ value_before_compactified_5
@test decompactifier(compactifier(value_before_compactified_6)) ≈ value_before_compactified_6
@test decompactifier(compactifier(value_before_compactified_7)) ≈ value_before_compactified_7

# setup, logistic
compactifier = Compactifier(lb, ub, ContinuousTransformations.Logistic())
decompactifier = Decompactifier(compactifier)

# values for try
value_before_compactified_1 = lb - 12.1     # out of bounds
value_before_compactified_2 = ub + 19.2     # out of bounds 
value_before_compactified_3 = lb            # lb
value_before_compactified_4 = ub            # ub
value_before_compactified_5 = (lb + ub) / 2 # within bounds
value_before_compactified_6 = lb - 4.1     # out of bounds, less extreme values
value_before_compactified_7 = ub + 1.2     # out of bounds, less extreme values

@test compactifier(-999999999.0) ≈ lb
@test compactifier(999999999.0) ≈ ub
@test decompactifier(lb) ≈ -Inf
@test decompactifier(ub) ≈ Inf
@test compactifier(0.0) ≈ (lb + ub) / 2

@test decompactifier(compactifier(value_before_compactified_1)) ≈ value_before_compactified_1
@test decompactifier(compactifier(value_before_compactified_2)) ≈ value_before_compactified_2
@test decompactifier(compactifier(value_before_compactified_3)) ≈ value_before_compactified_3
@test decompactifier(compactifier(value_before_compactified_4)) ≈ value_before_compactified_4
@test decompactifier(compactifier(value_before_compactified_5)) ≈ value_before_compactified_5
@test decompactifier(compactifier(value_before_compactified_6)) ≈ value_before_compactified_6
@test decompactifier(compactifier(value_before_compactified_7)) ≈ value_before_compactified_7
