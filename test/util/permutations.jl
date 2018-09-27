ps = get_all_permutations([[1.0]], 1)
@test [1.0] ∈ ps
@test length(ps) == 1

ps = get_all_permutations([[1.0], [3.0; 4.0]], 2)
@test [1.0; 3.0] ∈ ps
@test [1.0; 4.0] ∈ ps
@test length(ps) == 2

ps = get_all_permutations([[1.0; 2.0], [3.0]], 2)
@test [1.0; 3.0] ∈ ps
@test [2.0; 3.0] ∈ ps
@test length(ps) == 2

ps = get_all_permutations([[1.0; 2.0], [3.0; 4.0]], 2)
@test [1.0; 3.0] ∈ ps
@test [1.0; 4.0] ∈ ps
@test [2.0; 3.0] ∈ ps
@test [2.0; 4.0] ∈ ps
@test length(ps) == 4

ps = get_all_permutations([[1.0; 2.0], [3.0; 4.0; 5.0]], 2)
@test [1.0; 3.0] ∈ ps
@test [1.0; 4.0] ∈ ps
@test [1.0; 5.0] ∈ ps
@test [2.0; 3.0] ∈ ps
@test [2.0; 4.0] ∈ ps
@test [2.0; 5.0] ∈ ps
@test length(ps) == 6

ps = get_all_permutations([[1.0; 2.0], [3.0; 4.0], [5.0]], 3)
@test [1.0; 3.0; 5.0] ∈ ps
@test [1.0; 4.0; 5.0] ∈ ps
@test [2.0; 3.0; 5.0] ∈ ps
@test [2.0; 4.0; 5.0] ∈ ps
@test length(ps) == 4

ps = get_all_permutations([[1.0; 2.0], [3.0; 4.0], [5.0; 7.0; 9.0]], 3)
@test [1.0; 3.0; 5.0] ∈ ps
@test [1.0; 4.0; 5.0] ∈ ps
@test [2.0; 3.0; 5.0] ∈ ps
@test [2.0; 4.0; 5.0] ∈ ps
@test [1.0; 3.0; 7.0] ∈ ps
@test [1.0; 4.0; 7.0] ∈ ps
@test [2.0; 3.0; 7.0] ∈ ps
@test [2.0; 4.0; 7.0] ∈ ps
@test [1.0; 3.0; 9.0] ∈ ps
@test [1.0; 4.0; 9.0] ∈ ps
@test [2.0; 3.0; 9.0] ∈ ps
@test [2.0; 4.0; 9.0] ∈ ps
@test length(ps) == 12

ps = get_all_permutations([1.0:2.0, 3.0:4.0, 5.0:2.0:9.0], 3)
@test [1.0; 3.0; 5.0] ∈ ps
@test [1.0; 4.0; 5.0] ∈ ps
@test [2.0; 3.0; 5.0] ∈ ps
@test [2.0; 4.0; 5.0] ∈ ps
@test [1.0; 3.0; 7.0] ∈ ps
@test [1.0; 4.0; 7.0] ∈ ps
@test [2.0; 3.0; 7.0] ∈ ps
@test [2.0; 4.0; 7.0] ∈ ps
@test [1.0; 3.0; 9.0] ∈ ps
@test [1.0; 4.0; 9.0] ∈ ps
@test [2.0; 3.0; 9.0] ∈ ps
@test [2.0; 4.0; 9.0] ∈ ps
@test length(ps) == 12