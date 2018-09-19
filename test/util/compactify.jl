# setup
lb = 0.0
ub = 2.0
compactifier = Compactifier(lb, ub)
decompactifier = Compactifier(lb, ub)

# values for try
value_before_compactified_1 = lb - 12.1     # out of bounds
value_before_compactified_2 = ub + 19.2     # out of bounds 
value_before_compactified_3 = lb            # lb
value_before_compactified_4 = ub            # ub
value_before_compactified_5 = (lb + ub) / 2 # within bounds
value_before_compactified_6 = lb - 124.1     # out of bounds, more extreme values
value_before_compactified_7 = ub + 195.2     # out of bounds, more extreme values

@test compactifier(-999999999.0) ≈ lb
@test compactifier(999999999.0) ≈ ub
@test decompactify_approximately(compactifier, lb) ≈ -Inf
@test decompactify_approximately(compactifier, ub) ≈ Inf
@test compactifier(0.0) ≈ (lb + ub) / 2

@test decompactify_approximately(compactifier, compactifier(value_before_compactified_1)) ≈ value_before_compactified_1
@test decompactify_approximately(compactifier, compactifier(value_before_compactified_2)) ≈ value_before_compactified_2
@test decompactify_approximately(compactifier, compactifier(value_before_compactified_3)) ≈ value_before_compactified_3
@test decompactify_approximately(compactifier, compactifier(value_before_compactified_4)) ≈ value_before_compactified_4
@test decompactify_approximately(compactifier, compactifier(value_before_compactified_5)) ≈ value_before_compactified_5
@test decompactify_approximately(compactifier, compactifier(value_before_compactified_6)) ≈ value_before_compactified_6
@test decompactify_approximately(compactifier, compactifier(value_before_compactified_7)) ≈ value_before_compactified_7
