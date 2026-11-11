# Construction of modular subgroups from coset actions, their generators,
# membership tests and the S/T decomposition of SL(2,Z) matrices.
@testset "Basic modular group operations" begin
  s = @perm (2, 3)
  t = @perm (1, 2)
  G = modular_subgroup_via_right_action(s, t)
  GL = modular_subgroup_via_left_action(s, t)
  S = Oscar._perm_group(G)

  @test (@inferred index(G)) == 3
  @test (@inferred index(GL)) == 3

  @test (@inferred s_right_perm(GL)) == S(s^-1)
  @test (@inferred t_right_perm(GL)) == S(t^-1)

  s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
  t = cperm([1,4], [2,5,9,10,8], [3,7,6])
  G = modular_subgroup_via_right_action(s, t)
  @test (@inferred r_right_perm(G)) == cperm([1,7,9,10,6], [2,3], [4,5,8])
  @test (@inferred j_right_perm(G)) == cperm([1,8,3], [2,4,6], [5,7,10])

  # all words live in the same cached presentation
  SL2Z = @inferred s_t_presentation()
  @test s_t_presentation() === SL2Z
  @test parent(s_t_decomposition(matrix(ZZ, [2 1; 1 1]))) === SL2Z
  S, T = gens(SL2Z)
  actual_gens = @inferred word_gens(G)
  expected_gens = [
    S^-2,
    T^-2,
    S*T*S*T^-1*S^-1*T*S^-1,
    S*T^2*S*T^2*S^-1,
    S*T^5*S^-1
  ]
  @test all(is_word_element_of(w, G) for w in expected_gens)
  @test all(is_word_element_of(w, G) for w in actual_gens)

  actual_mat_gens = @inferred gens(G)
  expected_mat_gens = [[-1 0; 0 -1],
    [1 -2; 0 1],
    [2 -1; -3 2],
    [1 0; -5 1]
  ]
  @test all(m in G for m in actual_mat_gens)
  @test all(ZZMatrix(m) in G for m in expected_mat_gens)
  @test all(gen(G, i) in G for i in 1:length(actual_mat_gens))

  M = matrix(ZZ, [1 0; -4 1])
  @test (@inferred s_t_decomposition(M)) == S * T^4 * S^-1
  w = expected_gens[3]
  @test (@inferred Oscar._image_of_pt(w, G, 1)) == 1
  @test is_word_element_of(w, G)
end

@testset "Modular group equality, hashing and membership" begin
  s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
  t = cperm([1,4], [2,5,9,10,8], [3,7,6])
  G = modular_subgroup_via_right_action(s, t)

  # relabelling the cosets (fixing coset 1) yields the same subgroup
  c = cperm([2,3], [5,8,10])
  G2 = modular_subgroup_via_right_action(s^c, t^c)
  @test G == G2
  @test issubset(G, G2)
  @test hash(G) == hash(G2)
  @test length(Set([G, G2])) == 1

  # elements with equal matrices are equal and hash equally
  A = matrix(ZZ, [1 -2; 0 1])
  x = Oscar.ModularGroupElem(G, A)
  y = Oscar.ModularGroupElem(G, deepcopy(A))
  @test x == y
  @test hash(x) == hash(y)
  @test length(Set([x, y])) == 1
  @test x in G
  @test one(G) in G

  # elements can be created by calling the group, with a membership check by default
  @test (@inferred G(A)) == x
  @test parent(G(A)) === G
  @test G([1 -2; 0 1]) == x
  @test G([1, -2, 0, 1]) == x
  @test_throws ArgumentError G([1 1; 0 1])
  @test_throws ArgumentError G(identity_matrix(ZZ, 3))
  @test matrix(G([1 1; 0 1]; check = false)) == matrix(ZZ, [1 1; 0 1])
  H3 = modular_subgroup_via_right_action(cperm([2,3]), cperm([1,2]))
  @test parent(H3(G([1 0; 0 1]))) === H3
  @test_throws ArgumentError H3(G([1 1; 0 1]; check = false))

  # permutations from different symmetric groups are coerced
  H = modular_subgroup_via_right_action(cperm([2,3]), cperm([1,2]))
  @test index(H) == 3
  @test H != G
  S3 = symmetric_group(3)
  @test H == modular_subgroup_via_right_action(S3(cperm([2,3])), S3(cperm([1,2])))

  # the whole group SL(2,Z)
  W = modular_subgroup_via_right_action(cperm(), cperm())
  @test index(W) == 1
  @test matrix(ZZ, [1 0; 1 1]) in W
  @test isone(coset_right_action_of(A, W))

  # only 2x2 matrices of determinant 1 can be elements
  @test !(identity_matrix(ZZ, 3) in G)
  @test !(matrix(ZZ, [1 0; 0 -1]) in G)
  @test_throws ArgumentError s_t_decomposition(identity_matrix(ZZ, 3))
  @test_throws ArgumentError s_t_decomposition(matrix(ZZ, [1 0; 0 -1]))

  # membership must not be linear in the size of the matrix entries;
  # T^N lies in G iff N is a multiple of the length of the t-cycle of 1
  N = ZZ(2)^100
  @test matrix(ZZ, [1 N; 0 1]) in G
  @test !(matrix(ZZ, [1 N + 1; 0 1]) in G)

  # words are turned into matrices and coset actions correctly, also for huge exponents
  S, T = gens(s_t_presentation())
  @test all(Oscar._matrix_of_word(S^k) == matrix(ZZ, [0 -1; 1 0])^mod(k, 4) for k in -5:5)
  @test Oscar._matrix_of_word(T^N) == matrix(ZZ, [1 N; 0 1])
  M = matrix(ZZ, [N+1 N; 1 1])
  B = matrix(ZZ, [2 1; 1 1])
  @test Oscar._matrix_of_word(s_t_decomposition(M)) == M
  @test coset_right_action_of(matrix(ZZ, [1 N; 0 1]), G) == t_right_perm(G)^N
  @test coset_right_action_of(M * B, G) == coset_right_action_of(M, G) * coset_right_action_of(B, G)
  @test coset_left_action_of(M, G) == inv(coset_right_action_of(M, G))

  # order for some simple elements of SL(2,Z) must be correct
  x = Oscar.ModularGroupElem(W, matrix(ZZ, [1 0; 0 1]))
  @test order(x) == 1
  x = Oscar.ModularGroupElem(W, matrix(ZZ, [-1 0; 0 -1]))
  @test order(x) == 2
  x = Oscar.ModularGroupElem(W, matrix(ZZ, [0 -1; 1 0]))
  @test order(x) == 4

  # a nontrivial matrix with trace -2 has infinite order
  x = Oscar.ModularGroupElem(W, matrix(ZZ, [-1 1; 0 -1]))
  @test !is_finite_order(x)
  @test_throws Oscar.InfiniteOrderError order(x)

  # order must be infinite whenever the trace exceeds 2
  x = Oscar.ModularGroupElem(W, matrix(ZZ, [N -1; 1 0]))
  @test !is_finite_order(x)
  @test_throws Oscar.InfiniteOrderError order(x)
  # whenever the trace is 1, the order must be 6
  y = Oscar.ModularGroupElem(W, matrix(ZZ, [2 -3; 1 -1]))
  @test order(y) == 6

  @test_throws ArgumentError modular_subgroup_via_right_action(cperm([1,2,3]), cperm([1,2]))
  @test_throws ArgumentError modular_subgroup_via_left_action(cperm([1,2]), cperm([3,4]))

  # Powers of T belong to this subgroup exactly when the exponent is divisible by 3.
  @testset "Membership of powers of T" begin
    P = symmetric_group(3)
    G = modular_subgroup_via_right_action(one(P), P(cperm([1, 2, 3])))
    S, T = gens(s_t_presentation())

    for (n, expected) in [(-3, true), (-1, false), (0, true),
                          (1, false), (3, true)]
      A = matrix(ZZ, [1 n; 0 1])
      @test (A in G) == expected
      @test is_word_element_of(T^n, G) == expected
    end
  end

  # Left and right actions differ by inversion.
  @testset "Left and right coset actions" begin
    P = symmetric_group(3)
    s = one(P)
    t = P(cperm([1, 2, 3]))
    G = modular_subgroup_via_right_action(s, t)
    T = matrix(ZZ, [1 1; 0 1])

    @test index(G) == 3
    @test coset_right_action_of(T, G) == t
    @test coset_left_action_of(T, G) == t^-1

    H = modular_subgroup_via_left_action(s, t^-1)
    @test t_right_perm(H) == t
    @test H == G
  end

  # The i-th right coset representative maps the identity coset to coset i,
  # and the left coset representatives are their inverses.
  @testset "Coset representatives" begin
    s = cperm([1, 2], [3, 4], [5, 6], [7, 8], [9, 10])
    t = cperm([1, 4], [2, 5, 9, 10, 8], [3, 7, 6])
    G = modular_subgroup_via_right_action(s, t)
    S, T = gens(s_t_presentation())

    reps = @inferred right_coset_representatives(G)
    @test length(reps) == index(G)
    @test all(Oscar._image_of_pt(reps[i], G, 1) == i for i in 1:index(G))
    @test reps == [one(S), S, T*S, T, S*T, S*T*S, S*T^-1*S, S*T^-1, S*T^2, S*T^-2]

    @test right_coset_representatives(W) == [one(S)]

    lreps = @inferred left_coset_representatives(G)
    @test lreps == inv.(reps)
    @test all(Oscar._image_of_pt(inv(lreps[i]), G, 1) == i for i in 1:index(G))
    @test left_coset_representatives(W) == [one(S)]
  end

  # Equal indices do not imply equal subgroups or inclusion.
  @testset "Different subgroups of index three" begin
    P = symmetric_group(3)
    s = P(cperm([2, 3]))
    t = P(cperm([1, 2]))
    H = modular_subgroup_via_right_action(s, t)
    K = modular_subgroup_via_right_action(t, s)
    S = matrix(ZZ, [0 -1; 1 0])
    T = matrix(ZZ, [1 1; 0 1])

    @test index(H) == index(K) == 3
    @test S in H
    @test !(S in K)
    @test T in K
    @test !(T in H)
    @test H != K
    @test !issubset(H, K)
    @test !issubset(K, H)
  end

  # Relations and transitivity are separate requirements for a coset action.
  @testset "Invalid coset actions" begin
    P2 = symmetric_group(2)
    s = one(P2)
    t = P2(cperm([1, 2]))

    @test !defines_coset_action_s_t(s, t)
    @test_throws ArgumentError modular_subgroup_via_right_action(s, t)

    P3 = symmetric_group(3)
    u = P3(cperm([2, 3]))

    @test !defines_coset_action_s_t(u, u)
    @test_throws ArgumentError modular_subgroup_via_right_action(u, u)
  end

  # Invalid shapes and determinants cannot represent elements of SL(2,Z).
  @testset "Invalid matrices" begin
    bad_matrices = [
      identity_matrix(ZZ, 1),
      identity_matrix(ZZ, 3),
      matrix(ZZ, [1 0 0; 0 1 0]),
      zero_matrix(ZZ, 2, 2),
      matrix(ZZ, [1 0; 0 -1]),
      matrix(ZZ, [2 0; 0 1]),
    ]

    for A in bad_matrices
      @test !(A in W)
      @test !(A in H)
      @test_throws ArgumentError s_t_decomposition(A)
      @test_throws ArgumentError coset_right_action_of(A, H)
      @test_throws ArgumentError coset_left_action_of(A, H)
    end
  end

  # Finite-index modular subgroups are infinite, nonabelian and finitely generated.
  @testset "Group properties" begin
    for G in (W, H)
      @test !is_trivial(G)
      @test !is_finite(G)
      @test !is_abelian(G)
      @test !is_cyclic(G)
      @test is_finitely_generated(G)
    end
  end
end
