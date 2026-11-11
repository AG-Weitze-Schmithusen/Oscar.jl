include("exports.jl")
include("types.jl")

# Degree of the coset action, i.e. the index of the subgroup in SL(2,Z).
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
of ``{\rm SL}_2(\mathbb{Z})`` described by the permutations ``s`` and ``t``.

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

parent(x::ModularGroupElem) = x.parent
matrix(x::ModularGroupElem) = x.mat

one(G::ModularGroup) = ModularGroupElem(G, identity_matrix(ZZ, 2))
one(x::ModularGroupElem) = one(parent(x))

@doc raw"""
    (G::ModularGroup)(A::ZZMatrix; check::Bool=true)
    (G::ModularGroup)(A::AbstractVecOrMat; check::Bool=true)
    (G::ModularGroup)(x::ModularGroupElem; check::Bool=true)

Return the element of `G` given by the matrix `A` in ``{\rm SL}_2(\mathbb{Z})``,
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


function is_finite_order(x::ModularGroupElem)
  A = matrix(x)
  t = tr(A)
  return abs(t) < 2 || (t == 2 && isone(A)) || (t == -2 && isone(-A))
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
transitive coset action for ``{\rm SL}_2(\mathbb{Z})``.

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

@attr Int function index(G::ModularGroup)
  return _index(s_right_perm(G), t_right_perm(G))
end

const _SL2Z_FP_CACHE = Ref{FPGroup}()
const _MATRIX_HOM_CACHE = Ref{GAPGroupHomomorphism{FPGroup, MatGroup{ZZRingElem, ZZMatrix}}}()

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
``\langle S, T \mid S^4, (S^3 T)^3, S^2 T S^{-2} T^{-1} \rangle \cong {\rm SL}_2(\mathbb{Z})``,
where `S` and `T` correspond to the matrices
```math
S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \qquad T = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}.
```
All words in `S` and `T` returned by this package, for example by
[`s_t_decomposition`](@ref) and [`word_gens`](@ref), are elements of this group,
and words passed to [`is_word_element_of`](@ref) must be elements of it as well.
Every call returns the same group.

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

_matrix_S() = matrix(ZZ, [0 -1; 1 0])
_matrix_T() = matrix(ZZ, [1 1; 0 1])

# cache the homomorphism from the fp presentation of SL2Z to the matrix group
function _matrix_hom()
  if !isassigned(_MATRIX_HOM_CACHE)
    M = matrix_group([_matrix_S(), _matrix_T()])
    MS, MT = gens(M)

    SL2Z = _SL2Z_fp()
    _MATRIX_HOM_CACHE[] = hom(SL2Z, M, [MS, MT])
  end

  return _MATRIX_HOM_CACHE[]::GAPGroupHomomorphism{FPGroup, MatGroup{ZZRingElem, ZZMatrix}}
end

# cache the permutation group generated by s and t
@attr PermGroup function _perm_group(G::ModularGroup)
  return permutation_group(index(G), [s_right_perm(G), t_right_perm(G)])
end

# cache the homomorphism from SL2Z to P
@attr GAPGroupHomomorphism{FPGroup, PermGroup} function _coset_action_hom(G::ModularGroup)
  SL2Z = _SL2Z_fp()
  P = _perm_group(G)
  return hom(SL2Z, P, [s_right_perm(G), t_right_perm(G)])
end

inv(x::ModularGroupElem) = ModularGroupElem(parent(x), inv(matrix(x)))

function *(x::ModularGroupElem, y::ModularGroupElem)
  check_parent(x, y)
  return ModularGroupElem(parent(x), matrix(x) * matrix(y))
end

^(x::ModularGroupElem, n::Integer) =
  ModularGroupElem(parent(x), matrix(x)^n)

==(x::ModularGroupElem, y::ModularGroupElem) =
  matrix(x) == matrix(y)

Base.hash(x::ModularGroupElem, h::UInt) = hash(matrix(x), h)

@doc raw"""
    word_gens(G::ModularGroup)

Return generators of `G` as words in the generators `S` and `T` of
``{\rm SL}_2(\mathbb{Z})``, see [`s_t_decomposition`](@ref).

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
    right_coset_representatives(G::ModularGroup)

Return representatives of the right cosets of `G` in ``{\rm SL}_2(\mathbb{Z})``
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

Return representatives of the left cosets of `G` in ``{\rm SL}_2(\mathbb{Z})``
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
  w_gens = word_gens(G)
  phi = _matrix_hom()
  return [ModularGroupElem(G, matrix(phi(w))) for w in w_gens]
end

function gen(G::ModularGroup, i::Int)
  gs = gens(G)
  @req abs(i) <= length(gs) "i must be in the range -$(length(gs)):$(length(gs))"
  i > 0 && return gs[i]
  i < 0 && return inv(gs[-i])
  return one(G)
end

function number_of_generators(G::ModularGroup)
  return length(gens(G))
end

@doc raw"""
    s_t_decomposition(M::ZZMatrix)

Return the matrix `M` in ``{\rm SL}_2(\mathbb{Z})`` as a word in the generators
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

  MatS = _matrix_S()
  SL2Z = _SL2Z_fp()
  S, T = gens(SL2Z)

  decomp = one(SL2Z)

  # Euclidean algorithm on the second row: M * T^-k * S has second row
  # (d - k*c, -c), so the lower left entry decreases in absolute value.
  # T^{-k} = [1 -k; 0 1]. Do not power the matrix: fmpz_mat_pow takes a UInt64.
  while !iszero(M[2, 1])
    k = div(M[2, 2], M[2, 1])
    decomp = S^-1 * T^k * decomp
    M = M * matrix(ZZ, 2, 2, [1, -k, 0, 1]) * MatS
  end

  # now M[2, 1] = 0 and since det(M) == 1, we have M == +-T^r where r = +-M[1, 2]
  if isone(M[1, 1])
    decomp = T^M[1, 2] * decomp
  else
    decomp = S^2 * T^-M[1, 2] * decomp
  end

  return decomp
end

function coset_action_of(A::ZZMatrix, G::ModularGroup)
  w = s_t_decomposition(A)
  phi = _coset_action_hom(G)
  return phi(w)
end

@doc raw"""
    coset_right_action_of(A::ZZMatrix, G::ModularGroup)

Return the permutation describing the action of the matrix `A` in ``{\rm SL}_2(\mathbb{Z})``
on the right cosets of `G` by right multiplication.
"""
function coset_right_action_of(A::ZZMatrix, G::ModularGroup)
  return coset_action_of(A, G)
end

@doc raw"""
    coset_left_action_of(A::ZZMatrix, G::ModularGroup)

Return the permutation describing the action of the matrix `A` in ``{\rm SL}_2(\mathbb{Z})``
on the left cosets of `G` by left multiplication, see [`modular_subgroup_via_left_action`](@ref).
"""
function coset_left_action_of(A::ZZMatrix, G::ModularGroup)
  return coset_action_of(A, G)^-1
end

function _image_of_pt(w::FPGroupElem, G::ModularGroup, pt::Integer)
  @req 1 <= pt <= index(G) "Non-valid point"

  imgs = (s_right_perm(G), t_right_perm(G))

  for (i, e) in syllables(w) # i = 1 or 2, e = exponent
    pt = pt^(imgs[i]^e)
  end
  return pt
end

#TODO: could be optimized by tracing the image of the point 1 directly in the Euclidean loop,
# instead of decomposing A and take the resulting word apart again.
function Base.in(A::ZZMatrix, G::ModularGroup)
  (size(A) == (2, 2) && isone(det(A))) || return false
  w = s_t_decomposition(A)
  return _image_of_pt(w, G, 1) == 1
end

function Base.in(x::ModularGroupElem, G::ModularGroup)
  return parent(x) === G || matrix(x) in G
end

@doc raw"""
    is_word_element_of(w::FPGroupElem, G::ModularGroup)

Return whether the word `w` in the generators `S` and `T` of
``{\rm SL}_2(\mathbb{Z})`` represents an element of `G`.
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
  _check_compatible(parent(w), _SL2Z_fp())
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
