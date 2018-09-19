# setup
lb = 0.0
ub = 2.0
compactifier = PerlaTonettiWaugh.get_compactifier(lb, ub)
decompactifier = PerlaTonettiWaugh.get_decompactifier(lb, ub)

# values for try
value_before_compactified_1 = lb - 12.1     # out of bounds
value_before_compactified_2 = ub + 19.2     # out of bounds 
value_before_compactified_3 = lb            # lb
value_before_compactified_4 = ub            # ub
value_before_compactified_5 = (lb + ub) / 2 # within bounds

@test compactifier(-Inf) ≈ lb
@test compactifier(Inf) ≈ ub
@test decompactifier(lb) ≈ -Inf
@test decompactifier(ub) ≈ Inf
@test compactifier(0.0) ≈ (lb + ub) / 2

@test decompactifier(compactifier(value_before_compactified_1)) ≈ value_before_compactified_1
@test decompactifier(compactifier(value_before_compactified_2)) ≈ value_before_compactified_2
@test decompactifier(compactifier(value_before_compactified_3)) ≈ value_before_compactified_3
@test decompactifier(compactifier(value_before_compactified_4)) ≈ value_before_compactified_4
@test decompactifier(compactifier(value_before_compactified_5)) ≈ value_before_compactified_5