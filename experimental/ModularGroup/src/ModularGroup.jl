include("exports.jl")
include("types.jl")

# Degree of the coset action, i.e. the index of the subgroup in SL(2,Z) or PSL(2,Z).
# The trivial action (index 1) moves no points.
function _index(s::PermGroupElem, t::PermGroupElem)
  return max(1, largest_moved_point(s), largest_moved_point(t))
end

function modular_subgroup(s::PermGroupElem, t::PermGroupElem; check::Bool)
  # Coerce both permutations into the same symmetric group, so that the caller
  # may pass permutations living a priori in different parent symmetric groups.
  Sym = symmetric_group(_index(s, t))
  s = Sym(s)
  t = Sym(t)
  return ModularGroup(s, t; check)
end

@doc raw"""
    modular_subgroup_via_right_action(s::PermGroupElem, t::PermGroupElem; check::Bool=true)

Construct a `ModularGroup` object corresponding to the finite-index subgroup
of ``\mathrm{SL}_2(\mathbb{Z})`` described by the permutations ``s`` and ``t``.

For `check = true`, this
constructor tests if the given permutations actually describe the (right) coset action of the matrices
```math
S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \qquad T= \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix},
```
by checking that they act transitively and satisfy the relations
``s^4 = (s^3 t)^3 = s^2 t s^{-2} t^{-1} = 1``.
The point ``1`` corresponds to the coset of the identity matrix.

# Examples
```jldoctest
julia> s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
(1,2)(3,4)(5,6)(7,8)(9,10)

julia> t = cperm([1,4], [2,5,9,10,8], [3,7,6])
(1,4)(2,5,9,10,8)(3,7,6)

julia> G = modular_subgroup_via_right_action(s, t)
modular subgroup of index 10
```
"""
function modular_subgroup_via_right_action(s::PermGroupElem, t::PermGroupElem; check::Bool=true)
  return modular_subgroup(s, t; check)
end

@doc raw"""
    modular_subgroup_via_left_action(s::PermGroupElem, t::PermGroupElem; check::Bool=true)

Same as [`modular_subgroup_via_right_action`](@ref), but now the permutations describe the action
by left multiplication on the left cosets.

Under the bijection ``gH \mapsto Hg^{-1}`` this is the right action of the inverse matrices,
hence the subgroup is stored via the inverse permutations.

# Examples
```jldoctest
julia> s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
(1,2)(3,4)(5,6)(7,8)(9,10)

julia> t = cperm([1,4], [2,5,9,10,8], [3,7,6])
(1,4)(2,5,9,10,8)(3,7,6)

julia> G = modular_subgroup_via_left_action(s, t)
modular subgroup of index 10
```
"""
function modular_subgroup_via_left_action(s::PermGroupElem, t::PermGroupElem; check::Bool=true)
  return modular_subgroup(s^-1, t^-1; check)
end

function projective_modular_subgroup(s::PermGroupElem, t::PermGroupElem; check::Bool)
  Sym = symmetric_group(_index(s, t))
  s = Sym(s)
  t = Sym(t)
  return ProjectiveModularGroup(s, t; check)
end

@doc raw"""
    projective_modular_subgroup_via_right_action(s::PermGroupElem, t::PermGroupElem; check::Bool=true)

Construct a `ProjectiveModularGroup` object corresponding to the finite-index
subgroup of ``\mathrm{PSL}_2(\mathbb{Z})`` described by the permutations `s`
and `t`.

# Examples
```jldoctest
julia> projective_modular_subgroup_via_right_action(cperm([2, 3]), cperm([1, 2]))
projective modular subgroup of index 3
```
"""
function projective_modular_subgroup_via_right_action(s::PermGroupElem, t::PermGroupElem; check::Bool=true)
  Sym = symmetric_group(_index(s, t))
  s = Sym(s)
  t = Sym(t)
  return projective_modular_subgroup(s, t; check)
end

@doc raw"""
    projective_modular_subgroup_via_left_action(s::PermGroupElem, t::PermGroupElem; check::Bool=true)

Same as [`projective_modular_subgroup_via_right_action`](@ref), but now the permutations describe the action
by left multiplication on the left cosets.

Under the bijection ``gH \mapsto Hg^{-1}`` this is the right action of the inverse matrices,
hence the subgroup is stored via the inverse permutations.

# Examples
```jldoctest
julia> s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
(1,2)(3,4)(5,6)(7,8)(9,10)

julia> t = cperm([1,4], [2,5,9,10,8], [3,7,6])
(1,4)(2,5,9,10,8)(3,7,6)

julia> G = projective_modular_subgroup_via_left_action(s, t)
projective modular subgroup of index 10
```
"""
function projective_modular_subgroup_via_left_action(s::PermGroupElem, t::PermGroupElem; check::Bool=true)
  Sym = symmetric_group(_index(s, t))
  s = Sym(s)
  t = Sym(t)
  return projective_modular_subgroup(s^-1, t^-1; check)
end

function Base.:(==)(G::ModularGroup, H::ModularGroup)
  return index(G) == index(H) && issubset(G, H)
end

@attr NTuple{3, CycleType} function _cycle_structures(G::ModularGroup)
  return (cycle_structure(G.s), cycle_structure(G.t), cycle_structure(G.s * G.t))
end

function Base.hash(G::ModularGroup, h::UInt)
  s_struc, t_struc, s_t_struc = _cycle_structures(G)
  return hash((index(G), s_struc, t_struc, s_t_struc), h)
end

function Base.show(io::IO, G::ModularGroup)
  idx = index(G)
  print(io, "modular subgroup of index $(idx)")
end

function Base.show(io::IO, x::ModularGroupElem)
  A = matrix(x)
  idx = index(parent(x))
  print(io, "$(A) in modular group of index $(idx)")
end

elem_type(::Type{ModularGroup}) = ModularGroupElem
parent_type(::Type{ModularGroupElem}) = ModularGroup
elem_type(::Type{ProjectiveModularGroup}) = ProjectiveModularGroupElem
parent_type(::Type{ProjectiveModularGroupElem}) = ProjectiveModularGroup

parent(x::Union{ModularGroupElem, ProjectiveModularGroupElem}) = x.parent
matrix(x::Union{ModularGroupElem, ProjectiveModularGroupElem}) = x.mat

one(G::ModularGroup) = ModularGroupElem(G, identity_matrix(ZZ, 2))
one(G::ProjectiveModularGroup) = ProjectiveModularGroupElem(G, identity_matrix(ZZ, 2))
one(x::Union{ModularGroupElem, ProjectiveModularGroupElem}) = one(parent(x))

@doc raw"""
    (G::ModularGroup)(A::ZZMatrix; check::Bool=true)
    (G::ModularGroup)(A::AbstractVecOrMat; check::Bool=true)
    (G::ModularGroup)(x::ModularGroupElem; check::Bool=true)

Return the element of `G` given by the matrix `A` in ``\mathrm{SL}_2(\mathbb{Z})``,
or by the matrix of `x`, which may be an element of a different modular subgroup.

For `check = true`, an `ArgumentError` is thrown if the matrix is not in `G`.
For `check = false`, no test is performed.

# Examples
```jldoctest
julia> G = modular_subgroup_via_right_action(cperm([2, 3]), cperm([1, 2]))
modular subgroup of index 3

julia> x = G([1 2; 0 1])
[1 2; 0 1] in modular group of index 3

julia> x == G(matrix(ZZ, [1 2; 0 1]))
true

julia> G([1 1; 0 1])
ERROR: ArgumentError: Element not in the group
```
"""
function (G::ModularGroup)(A::ZZMatrix; check::Bool=true)
  if check
    @req A in G "Element not in the group"
  end
  return ModularGroupElem(G, A)
end

function (G::ModularGroup)(A::AbstractVecOrMat; check::Bool=true)
  return G(matrix(ZZ, 2, 2, A); check)
end

function (G::ModularGroup)(x::ModularGroupElem; check::Bool=true)
  return G(matrix(x); check)
end

@doc raw"""
    (G::ProjectiveModularGroup)(A::ZZMatrix; check::Bool=true)
    (G::ProjectiveModularGroup)(A::AbstractVecOrMat; check::Bool=true)
    (G::ProjectiveModularGroup)(x::ProjectiveModularGroupElem; check::Bool=true)

Return the element of `G` given by the matrix `A` in ``\mathrm{PSL}_2(\mathbb{Z})``,
or by the matrix of `x`, which may be an element of a different projective modular
subgroup.

Note that `±A` is the same element of ``\mathrm{PSL}_2(\mathbb{Z})``. Thus
`A = [a b; c d]` is normalized such that `c > 0`, or `c = 0` and `d > 0`.

For `check = true`, an `ArgumentError` is thrown if the matrix is not in `G`.
For `check = false`, no test is performed.

# Examples
```jldoctest
julia> G = projective_modular_subgroup_via_right_action(cperm([2, 3]), cperm([1, 2]))
projective modular subgroup of index 3

julia> x = G([1 2; 0 1])
[1 2; 0 1] in projective modular group of index 3

julia> y = G([-1 -2; 0 -1])
[1 2; 0 1] in projective modular group of index 3
```
"""
function (G::ProjectiveModularGroup)(A::ZZMatrix; check::Bool=true)
  if check
    @req A in G "Element not in the group"
  end
  return ProjectiveModularGroupElem(G, A)
end

function (G::ProjectiveModularGroup)(A::AbstractVecOrMat; check::Bool=true)
  return G(matrix(ZZ, 2, 2, A); check)
end

function (G::ProjectiveModularGroup)(x::ProjectiveModularGroupElem; check::Bool=true)
  return G(matrix(x); check)
end

function order(::Type{T}, x::ModularGroupElem) where {T<:IntegerUnion}
  A = matrix(x)
  t = tr(A)

  if t == 2
    isone(A) || throw(InfiniteOrderError(x))
    return T(1)
  elseif t == -2
    isone(-A) || throw(InfiniteOrderError(x))
    return T(2)
  elseif t == -1
    return T(3)
  elseif t == 0
    return T(4)
  elseif t == 1
    return T(6)
  else
    throw(InfiniteOrderError(x))
  end
end

function order(::Type{T}, x::ProjectiveModularGroupElem) where {T<:IntegerUnion}
  A = matrix(x)
  t = tr(A)

  if t == 2
    isone(A) || throw(InfiniteOrderError(x))
    return T(1)
  elseif t == 0
    return T(2)
  elseif abs(t) == 1
    return T(3)
  else
    throw(InfiniteOrderError(x))
  end
end

function is_finite_order(x::ModularGroupElem)
  A = matrix(x)
  t = tr(A)
  return abs(t) < 2 || (t == 2 && isone(A)) || (t == -2 && isone(-A))
end

function is_finite_order(x::ProjectiveModularGroupElem)
  A = matrix(x)
  t = tr(A)
  return abs(t) < 2 || (t == 2 && isone(A))
end

@doc raw"""
    s_right_perm(G::ModularGroup)

Return the permutation describing the action of the matrix ``S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}``
on the right cosets of `G`.
"""
function s_right_perm(G::ModularGroup)
  return G.s
end

@doc raw"""
    t_right_perm(G::ModularGroup)

Return the permutation describing the action of the matrix ``T= \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}``
on the right cosets of `G`.
"""
function t_right_perm(G::ModularGroup)
  return G.t
end

@doc raw"""
    r_right_perm(G::ModularGroup)

Return the permutation describing the action of the matrix ``R= \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}``
on the right cosets of `G`.
"""
@attr PermGroupElem function r_right_perm(G::ModularGroup)
  return G.s^-1*G.t^-1*G.s
end

@doc raw"""
    j_right_perm(G::ModularGroup)

Return the permutation describing the action of the matrix ``J= \begin{pmatrix} 0 & 1 \\ -1 & 1 \end{pmatrix}``
on the right cosets of `G`.
"""
@attr PermGroupElem function j_right_perm(G::ModularGroup)
  return G.s^-1*G.t^-1
end

@doc raw"""
    defines_coset_action_s_t(s::PermGroupElem, t::PermGroupElem)

Return whether the elements `s` and `t` of a permutation group define a valid 
transitive coset action for ``\mathrm{SL}_2(\mathbb{Z})``.

# Examples
```jldoctest
julia> Sym = symmetric_group(3);

julia> s = Sym(cperm([2, 3])); t = Sym(cperm([1, 2]));

julia> defines_coset_action_s_t(s, t)
true

julia> s_bad = Sym(cperm([1, 2, 3])); # Order 3 contradicts S^4 = 1

julia> defines_coset_action_s_t(s_bad, t)
false
```
"""
function defines_coset_action_s_t(s::PermGroupElem, t::PermGroupElem)
  isone(s^4) || return false
  isone((s^3*t)^3) || return false
  isone(s^2*t*s^-2*t^-1) || return false
  return is_transitive(permutation_group(_index(s, t), [s, t]))
end

@doc raw"""
    defines_projective_coset_action(s::PermGroupElem, t::PermGroupElem)

Return whether the elements `s` and `t` of a permutation group define a valid
transitive coset action for ``{\rm PSL}_2(\mathbb{Z})``.

# Examples
```jldoctest
julia> s = cperm([2, 3]); t = cperm([1, 2]);

julia> defines_projective_coset_action(s, t)
true

julia> s_bad = cperm([1, 2, 3]); # Order 3 contradicts S^2 = 1

julia> defines_projective_coset_action(s_bad, t)
false
```
"""
function defines_projective_coset_action(s::PermGroupElem, t::PermGroupElem)
  isone(s^2) || return false
  isone((s*t)^3) || return false

  return is_transitive(permutation_group(_index(s, t), [s, t]))
end

@attr Int function index(G::ModularGroup)
  return _index(s_right_perm(G), t_right_perm(G))
end

const _SL2Z_FP_CACHE = Ref{FPGroup}()

# cache the SL2Z presentation so that it can be reused consistently
# (especially in testing)
function _SL2Z_fp()
  if !isassigned(_SL2Z_FP_CACHE)
    # Words like T^k with huge k arise from matrices with large entries;
    # the syllable representation stores them in constant space.
    F = free_group([:S, :T]; eltype = :syllable)
    S, T = gens(F)

    SL2Z, _ = quo(F, [S^4, (S^3*T)^3, S^2*T*S^-2*T^-1])

    _SL2Z_FP_CACHE[] = SL2Z
  end

  return _SL2Z_FP_CACHE[]::FPGroup
end

@doc raw"""
    s_t_presentation()

Return the finitely presented group
``\langle S, T \mid S^4, (S^3 T)^3, S^2 T S^{-2} T^{-1} \rangle \cong \mathrm{SL}_2(\mathbb{Z})``,
where `S` and `T` correspond to the matrices
```math
S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \qquad T = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}.
```
All words in `S` and `T` returned by this package, for example by
[`s_t_decomposition`](@ref) and [`word_gens`](@ref), are elements of this group or of the group returned
by [`projective_s_t_presentation`](@ref), and words passed to [`is_word_element_of`](@ref) must be elements
of one of those groups as well.

# Examples
```jldoctest
julia> SL2Z = s_t_presentation()
Finitely presented group

julia> S, T = gens(SL2Z);

julia> s_t_decomposition(matrix(ZZ, [1 2; 0 1])) == T^2
true
```
"""
s_t_presentation() = _SL2Z_fp()

const _PSL2Z_FP_CACHE = Ref{FPGroup}()
const _SL2Z_TO_PSL2Z_CACHE = Ref{GAPGroupHomomorphism{FPGroup, FPGroup}}()

# cache the presentation <S, T | S^2, (S*T)^3> of PSL2Z, analogous to `_SL2Z_fp`
function _PSL2Z_fp()
  if !isassigned(_PSL2Z_FP_CACHE)
    F = free_group([:S, :T]; eltype = :syllable)
    S, T = gens(F)

    PSL2Z, _ = quo(F, [S^2, (S*T)^3])

    _PSL2Z_FP_CACHE[] = PSL2Z
  end

  return _PSL2Z_FP_CACHE[]::FPGroup
end

@doc raw"""
    projective_s_t_presentation()

Return the finitely presented group
``\langle S, T \mid S^2, (S T)^3 \rangle \cong \mathrm{PSL}_2(\mathbb{Z})``, where `S` and `T`
correspond to the matrices
```math
S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \qquad T = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}.
```
All words in `S` and `T` returned by this package, for example by
[`s_t_decomposition`](@ref) and [`word_gens`](@ref), are elements of this group or of the group returned
by [`s_t_presentation`](@ref), and words passed to [`is_word_element_of`](@ref) must be elements of one of
those groups as well.

# Examples
```jldoctest
julia> PSL2Z = projective_s_t_presentation()
Finitely presented group

julia> S, T = gens(PSL2Z);

julia> sl2z_to_psl2z()(s_t_decomposition(matrix(ZZ, [0 1; -1 0]))) == S
true
```
"""
projective_s_t_presentation() = _PSL2Z_fp()

# cache the projection from `_SL2Z_fp()` onto `_PSL2Z_fp()`, mapping S to S and T to T
function sl2z_to_psl2z()
  if !isassigned(_SL2Z_TO_PSL2Z_CACHE)
    # This is a homomorphism since S^2 = 1 turns every relator of SL2Z into a
    # consequence of S^2 and (S*T)^3. Checking it in GAP would require solving
    # the word problem in an infinite finitely presented group, so skip it.
    _SL2Z_TO_PSL2Z_CACHE[] = hom(_SL2Z_fp(), _PSL2Z_fp(), gens(_PSL2Z_fp()); check = false)
  end

  return _SL2Z_TO_PSL2Z_CACHE[]::GAPGroupHomomorphism{FPGroup, FPGroup}
end

# The following helpers evaluate a word in S and T (an element of `s_t_presentation()`)
# syllable by syllable, so a syllable T^k with huge k costs as much as T.

# Matrix of S^e (i = 1) or T^e (i = 2) in closed form. Powering a ZZMatrix
# directly fails for exponents that do not fit into a UInt.
function _matrix_of_syllable(i::Int, e::IntegerUnion)
  matrix_S = matrix(ZZ, [0 -1; 1 0])
  i == 2 && return matrix(ZZ, 2, 2, [1, e, 0, 1])
  r = mod(e, 4)
  r == 0 && return identity_matrix(ZZ, 2)
  r == 1 && return matrix_S
  r == 2 && return -identity_matrix(ZZ, 2)
  return - matrix_S # S^3 = S^-1 = -S
end

function _matrix_of_word(w::FPGroupElem)
  M = identity_matrix(ZZ, 2)
  for (i, e) in syllables(w)
    M = M * _matrix_of_syllable(i, e)
  end
  return M
end

# Permutation by which the word w acts on the right cosets of G.
function _perm_of_word(w::FPGroupElem, G::Union{ModularGroup, ProjectiveModularGroup})
  imgs = (s_right_perm(G), t_right_perm(G))
  p = one(imgs[1])
  for (i, e) in syllables(w)
    p = p * imgs[i]^e
  end
  return p
end

# cache the permutation group generated by s and t
@attr PermGroup function _perm_group(G::Union{ModularGroup, ProjectiveModularGroup})
  return permutation_group(index(G), [s_right_perm(G), t_right_perm(G)])
end

# cache the homomorphism from SL2Z to P; only needed for `preimage` in `word_gens`
@attr GAPGroupHomomorphism{FPGroup, PermGroup} function _coset_action_hom(G::ModularGroup)
  SL2Z = _SL2Z_fp()
  P = _perm_group(G)
  return hom(SL2Z, P, [s_right_perm(G), t_right_perm(G)])
end

# cache the homomorphism from PSL2Z to P; only needed for `preimage` in `word_gens`
@attr GAPGroupHomomorphism{FPGroup, PermGroup} function _coset_action_hom(G::ProjectiveModularGroup)
  PSL2Z = _PSL2Z_fp()
  P = _perm_group(G)
  return hom(PSL2Z, P, [s_right_perm(G), t_right_perm(G)])
end

inv(x::ModularGroupElem) = ModularGroupElem(parent(x), inv(matrix(x)))

inv(x::ProjectiveModularGroupElem) = ProjectiveModularGroupElem(parent(x), inv(matrix(x)))

function *(x::ModularGroupElem, y::ModularGroupElem)
  check_parent(x, y)
  return ModularGroupElem(parent(x), matrix(x) * matrix(y))
end

function *(x::ProjectiveModularGroupElem, y::ProjectiveModularGroupElem)
  check_parent(x, y)
  return ProjectiveModularGroupElem(parent(x), matrix(x) * matrix(y))
end

^(x::ModularGroupElem, n::Integer) =
  ModularGroupElem(parent(x), matrix(x)^n)

^(x::ProjectiveModularGroupElem, n::Integer) =
  ProjectiveModularGroupElem(parent(x), matrix(x)^n)

==(x::ModularGroupElem, y::ModularGroupElem) =
  matrix(x) == matrix(y)

==(x::ProjectiveModularGroupElem, y::ProjectiveModularGroupElem) =
  matrix(x) == matrix(y)

Base.hash(x::Union{ModularGroupElem, ProjectiveModularGroupElem}, h::UInt) = hash(matrix(x), h)

@doc raw"""
    word_gens(G::ModularGroup)

Return generators of `G` as words in the generators `S` and `T` of
``\mathrm{SL}_2(\mathbb{Z})``, see [`s_t_presentation`](@ref).

# Examples
```jldoctest
julia> G = modular_subgroup_via_right_action(cperm([2, 3]), cperm([1, 2]))
modular subgroup of index 3

julia> word_gens(G)
2-element Vector{FPGroupElem}:
 S
 T^-2
```
"""
function word_gens(G::ModularGroup)
  P = _perm_group(G)
  phi = _coset_action_hom(G)

  Hperm, _ = stabilizer(P, 1)

  H, inc = preimage(phi, Hperm)

  return [inc(h) for h in gens(H)]
end

@doc raw"""
    word_gens(G::ProjectiveModularGroup)

Return generators of `G` as words in the generators `S` and `T` of
``\mathrm{PSL}_2(\mathbb{Z})``, see [`projective_s_t_presentation`](@ref).

# Examples
```jldoctest
julia> G = projective_modular_subgroup_via_right_action(cperm([2, 3]), cperm([1, 2]))
projective modular subgroup of index 3

julia> word_gens(G)
2-element Vector{FPGroupElem}:
 S
 T^-2
```
"""
function word_gens(G::ProjectiveModularGroup)
  P = _perm_group(G)
  phi = _coset_action_hom(G)

  Hperm, _ = stabilizer(P, 1)

  H, inc = preimage(phi, Hperm)

  return [inc(h) for h in gens(H)]
end

@doc raw"""
    right_coset_representatives(G::ModularGroup)

Return representatives of the right cosets of `G` in ``\mathrm{SL}_2(\mathbb{Z})``
as words in the generators `S` and `T`, see [`s_t_decomposition`](@ref).

The ``i``-th entry represents the coset with label ``i``, i.e. it maps the coset
``1`` of the identity to the coset ``i`` under the action given by
[`s_right_perm`](@ref) and [`t_right_perm`](@ref).
The representatives are found by a breadth-first search through the coset graph,
so each one is a shortest word in `S`, `T` and their inverses.

# Examples
```jldoctest
julia> G = modular_subgroup_via_right_action(cperm([1, 2]), cperm([2, 3]))
modular subgroup of index 3

julia> right_coset_representatives(G)
3-element Vector{FPGroupElem}:
 <identity ...>
 S
 S*T
```
"""
@attr Vector{FPGroupElem} function right_coset_representatives(G::ModularGroup)
  n = index(G)
  S, T = gens(_SL2Z_fp())
  perms = (s_right_perm(G), s_right_perm(G)^-1, t_right_perm(G), t_right_perm(G)^-1)
  words = (S, S^-1, T, T^-1)

  reps = Vector{FPGroupElem}(undef, n)
  reps[1] = one(S)
  queue = [1]
  sizehint!(queue, n)
  head = 1

  while head <= length(queue)
    x = queue[head]
    head += 1
    for (p, w) in zip(perms, words)
      y = x^p
      isassigned(reps, y) && continue
      reps[y] = reps[x] * w
      push!(queue, y)
    end
  end

  return reps
end

@doc raw"""
    left_coset_representatives(G::ModularGroup)

Return representatives of the left cosets of `G` in ``\mathrm{SL}_2(\mathbb{Z})``
as words in the generators `S` and `T`, see [`s_t_decomposition`](@ref).

These are the inverses of the words returned by [`right_coset_representatives`](@ref),
so the ``i``-th entry represents the left coset with label ``i`` under the bijection
``gG \mapsto Gg^{-1}``, see [`coset_left_action_of`](@ref).

# Examples
```jldoctest
julia> G = modular_subgroup_via_right_action(cperm([1, 2]), cperm([2, 3]))
modular subgroup of index 3

julia> left_coset_representatives(G)
3-element Vector{FPGroupElem}:
 <identity ...>
 S^-1
 T^-1*S^-1
```
"""
@attr Vector{FPGroupElem} function left_coset_representatives(G::ModularGroup)
  return [inv(w) for w in right_coset_representatives(G)]
end

@attr Vector{ModularGroupElem} function gens(G::ModularGroup)
  return [ModularGroupElem(G, _matrix_of_word(w)) for w in word_gens(G)]
end

@attr Vector{ProjectiveModularGroupElem} function gens(G::ProjectiveModularGroup)
  return [ProjectiveModularGroupElem(G, _matrix_of_word(w)) for w in word_gens(G)]
end

function gen(G::Union{ModularGroup, ProjectiveModularGroup}, i::Int)
  gs = gens(G)
  @req abs(i) <= length(gs) "i must be in the range -$(length(gs)):$(length(gs))"
  i > 0 && return gs[i]
  i < 0 && return inv(gs[-i])
  return one(G)
end

function number_of_generators(G::Union{ModularGroup, ProjectiveModularGroup})
  return length(gens(G))
end

@doc raw"""
    s_t_decomposition(M::ZZMatrix)

Return the matrix `M` in ``\mathrm{SL}_2(\mathbb{Z})`` as a word in the generators
```math
S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \qquad T = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}.
```
The word is an element of the finitely presented group returned by
[`s_t_presentation`](@ref).

# Examples
```jldoctest
julia> G = modular_subgroup_via_right_action(cperm([2, 3]), cperm([1, 2]))
modular subgroup of index 3

julia> A = gens(G)[2]
[1 -2; 0 1] in modular group of index 3

julia> s_t_decomposition(matrix(A))
T^-2
```
"""
function s_t_decomposition(M::ZZMatrix)
  @req size(M) == (2, 2) "Matrix needs to be in SL(2, Z)"
  @req isone(det(M)) "Matrix needs to be in SL(2, Z)"

  SL2Z = _SL2Z_fp()
  S, T = gens(SL2Z)

  decomp = one(SL2Z)

  # Euclidean algorithm on the second row: M * T^-k * S has second row
  # (d - k*c, -c), so the lower left entry decreases in absolute value.
  # T^{-k} = [1 -k; 0 1]. Do not power the matrix: fmpz_mat_pow takes a UInt64.
  while !iszero(M[2, 1])
    k = div(M[2, 2], M[2, 1])
    decomp = S^-1 * T^k * decomp
    M = M * matrix(ZZ, 2, 2, [-k -1; 1 0])
  end

  # now M[2, 1] = 0 and since det(M) == 1, we have M == +-T^r where r = +-M[1, 2]
  if isone(M[1, 1])
    decomp = T^M[1, 2] * decomp
  else
    decomp = S^2 * T^-M[1, 2] * decomp
  end

  return decomp
end

function coset_action_of(A::ZZMatrix, G::Union{ModularGroup, ProjectiveModularGroup})
  return _perm_of_word(s_t_decomposition(A), G)
end

@doc raw"""
    coset_right_action_of(A::ZZMatrix, G::ModularGroup)

Return the permutation describing the action of the matrix `A` in ``\mathrm{SL}_2(\mathbb{Z})``
on the right cosets of `G` by right multiplication.
"""
function coset_right_action_of(A::ZZMatrix, G::ModularGroup)
  return coset_action_of(A, G)
end

@doc raw"""
    coset_right_action_of(A::ZZMatrix, G::ProjectiveModularGroup)

Return the permutation describing the action of the matrix `A` in ``\mathrm{PSL}_2(\mathbb{Z})``
on the right cosets of `G` by right multiplication.
"""
function coset_right_action_of(A::ZZMatrix, G::ProjectiveModularGroup)
  return coset_action_of(A, G)
end


@doc raw"""
    coset_left_action_of(A::ZZMatrix, G::ModularGroup)

Return the permutation describing the action of the matrix `A` in ``\mathrm{SL}_2(\mathbb{Z})``
on the left cosets of `G` by left multiplication, see [`modular_subgroup_via_left_action`](@ref).
"""
function coset_left_action_of(A::ZZMatrix, G::ModularGroup)
  return coset_action_of(A, G)^-1
end

@doc raw"""
    coset_left_action_of(A::ZZMatrix, G::ProjectiveModularGroup)

Return the permutation describing the action of the matrix `A` in ``\mathrm{PSL}_2(\mathbb{Z})``
on the left cosets of `G` by left multiplication, see [`projective_modular_subgroup_via_left_action`](@ref).
"""
function coset_left_action_of(A::ZZMatrix, G::ProjectiveModularGroup)
  return coset_action_of(A, G)^-1
end

function _image_of_pt(w::FPGroupElem, G::Union{ModularGroup, ProjectiveModularGroup}, pt::Integer)
  @req 1 <= pt <= index(G) "Non-valid point"

  imgs = (s_right_perm(G), t_right_perm(G))

  for (i, e) in syllables(w) # i = 1 or 2, e = exponent
    pt = pt^(imgs[i]^e)
  end
  return pt
end

#TODO: could be optimized by tracing the image of the point 1 directly in the Euclidean loop,
# instead of decomposing A and take the resulting word apart again.
function Base.in(A::ZZMatrix, G::Union{ModularGroup, ProjectiveModularGroup})
  (size(A) == (2, 2) && isone(det(A))) || return false
  w = s_t_decomposition(A)
  return _image_of_pt(w, G, 1) == 1
end

function Base.in(x::ModularGroupElem, G::ModularGroup)
  return parent(x) === G || matrix(x) in G
end

function Base.in(x::ProjectiveModularGroupElem, G::ProjectiveModularGroup)
  return parent(x) === G || matrix(x) in G
end


@doc raw"""
    is_word_element_of(w::FPGroupElem, G::ModularGroup)

Return whether the word `w` in the generators `S` and `T` represents an element of `G`.

The word must be an element of the group returned by [`s_t_presentation`](@ref).

# Examples
```jldoctest
julia> G = modular_subgroup_via_right_action(cperm([2, 3]), cperm([1, 2]))
modular subgroup of index 3

julia> S, T = gens(s_t_presentation());

julia> is_word_element_of(T^2, G)
true

julia> is_word_element_of(T, G)
false
```
"""
function is_word_element_of(w::FPGroupElem, G::ModularGroup)
  P = parent(w)
  @req _check_compatible(P, _SL2Z_fp()) "w must be a word in the presentation of SL(2,Z)"
  return _image_of_pt(w, G, 1) == 1
end

@doc raw"""
    is_word_element_of(w::FPGroupElem, G::ProjectiveModularGroup)

Return whether the word `w` in the generators `S` and `T` represents an element of `G`.

The word must be an element of the group returned by [`s_t_presentation`](@ref) or by
[`projective_s_t_presentation`](@ref). A word in the latter determines a matrix only up to
sign, so it is accepted only if ``-I`` acts trivially on the cosets of `G`, that is, if `G`
is a projective modular subgroup or a modular subgroup containing ``-I``.

# Examples
```jldoctest
julia> G = projective_modular_subgroup_via_right_action(cperm([2, 3]), cperm([1, 2]))
projective modular subgroup of index 3

julia> S, T = gens(projective_s_t_presentation());

julia> is_word_element_of(S, G)
true
```
"""
function is_word_element_of(w::FPGroupElem, G::ProjectiveModularGroup)
  P = parent(w)
  if !_check_compatible(P, _SL2Z_fp(); error = false)
    @req _check_compatible(P, _PSL2Z_fp(); error = false) "w must be a word in the presentation of SL(2,Z) or PSL(2,Z)"
    # s^2 is the action of -I. If it is trivial, s and t satisfy S^2 = (S*T)^3 = 1,
    # so they define an action of PSL(2,Z) and words in its presentation can be traced.
    @req isone(s_right_perm(G)^2) "words in PSL(2,Z) need a subgroup containing -I"
  end
  return _image_of_pt(w, G, 1) == 1
end

function Base.issubset(H::ModularGroup, G::ModularGroup)
  nH = index(H)
  nG = index(G)

  # Necessary numerical conditions: [SL2Z : H] = [SL2Z : G] * [G : H]
  nH >= nG || return false
  iszero(mod(nH, nG)) || return false

  sH = s_right_perm(H)
  tH = t_right_perm(H)
  sG = s_right_perm(G)
  tG = t_right_perm(G)

  f = zeros(Int, nH) # f[x] = image of H-coset x in the G-cosets
  f[1] = 1 # identity coset |-> identity coset

  queue = Vector{Int}(undef, nH)
  queue[1] = 1
  head = 1
  tail = 1

  @inbounds while head <= tail
    x  = queue[head]
    head += 1
    fx = f[x]

    for (pH, pG) in ((sH, sG), (tH, tG))
      y  = x^pH  # neighbour in the H-coset graph
      fy = fx^pG # forced image under equivariance

      if fy < 1 || fy > nG
        return false
      elseif f[y] == 0
        f[y] = fy
        tail += 1
        queue[tail] = y
      elseif f[y] != fy
        return false # conflicting assignment
      end
    end
  end

  # Transitivity of the H-action (already guaranteed by the constructor,
  # but cheap to assert that the map was defined everywhere).
  return tail == nH
end

is_abelian(::ModularGroup) = false
is_finite(::ModularGroup) = false
is_cyclic(::ModularGroup) = false
is_trivial(::ModularGroup) = false
is_finitely_generated(::ModularGroup) = true
#TODO random elements?

# Projective modular subgroups are represented by transitive permutation
# actions of the presentation of PSL(2, Z).
function Base.:(==)(G::ProjectiveModularGroup, H::ProjectiveModularGroup)
  return index(G) == index(H) && issubset(G, H)
end

@attr NTuple{3, CycleType} function _cycle_structures(G::ProjectiveModularGroup)
  return (cycle_structure(G.s), cycle_structure(G.t), cycle_structure(G.s * G.t))
end

function Base.hash(G::ProjectiveModularGroup, h::UInt)
  s_struc, t_struc, s_t_struc = _cycle_structures(G)
  return hash((index(G), s_struc, t_struc, s_t_struc), h)
end

function Base.show(io::IO, G::ProjectiveModularGroup)
  idx = index(G)
  print(io, "projective modular subgroup of index $(idx)")
end

function Base.show(io::IO, x::ProjectiveModularGroupElem)
  A = matrix(x)
  idx = index(parent(x))
  print(io, "$(A) in projective modular group of index $(idx)")
end

@doc raw"""
    s_right_perm(G::ProjectiveModularGroup)

Return the permutation describing the action of the matrix ``S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}``
on the right cosets of `G`.
"""
function s_right_perm(G::ProjectiveModularGroup)
  return G.s
end

@doc raw"""
    t_right_perm(G::ProjectiveModularGroup)

Return the permutation describing the action of the matrix ``T= \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}``
on the right cosets of `G`.
"""
function t_right_perm(G::ProjectiveModularGroup)
  return G.t
end

@doc raw"""
    r_right_perm(G::ProjectiveModularGroup)

Return the permutation describing the action of the matrix ``R= \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}``
on the right cosets of `G`.
"""
@attr PermGroupElem function r_right_perm(G::ProjectiveModularGroup)
  return G.s^-1*G.t^-1*G.s
end

@doc raw"""
    j_right_perm(G::ProjectiveModularGroup)

Return the permutation describing the action of the matrix ``J= \begin{pmatrix} 0 & 1 \\ -1 & 1 \end{pmatrix}``
on the right cosets of `G`.
"""
@attr PermGroupElem function j_right_perm(G::ProjectiveModularGroup)
  return G.s^-1*G.t^-1
end

@attr Int function index(G::ProjectiveModularGroup)
  return _index(s_right_perm(G), t_right_perm(G))
end

function Base.issubset(H::ProjectiveModularGroup, G::ProjectiveModularGroup)
  nH = index(H)
  nG = index(G)

  # Necessary numerical conditions
  nH >= nG || return false
  iszero(mod(nH, nG)) || return false

  sH = s_right_perm(H)
  tH = t_right_perm(H)
  sG = s_right_perm(G)
  tG = t_right_perm(G)

  f = zeros(Int, nH)
  f[1] = 1

  queue = Vector{Int}(undef, nH)
  queue[1] = 1
  head = 1
  tail = 1

  @inbounds while head <= tail
    x = queue[head]
    head += 1
    fx = f[x]

    for (pH, pG) in ((sH, sG), (tH, tG))
      y = x^pH
      fy = fx^pG

      if fy < 1 || fy > nG
        return false
      elseif f[y] == 0
        f[y] = fy
        tail += 1
        queue[tail] = y
      elseif f[y] != fy
        return false
      end
    end
  end

  return tail == nH
end

is_abelian(::ProjectiveModularGroup) = false
is_finite(::ProjectiveModularGroup) = false
is_cyclic(::ProjectiveModularGroup) = false
is_trivial(::ProjectiveModularGroup) = false
is_finitely_generated(::ProjectiveModularGroup) = true
