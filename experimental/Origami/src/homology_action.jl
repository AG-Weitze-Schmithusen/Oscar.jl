@doc raw"""
    homology(o::Origami)

Return a ``\mathbb{Z}``-basis of the homology ``H = H_1(X, \mathbb{Z})`` of the origami `o`.

The homology is computed as cellular homology with the ``d`` squares of `o` as 2-cells. Let ``\sigma_i`` be the bottom
edge and ``\xi_i`` the left edge of the ``i``-th square. An element of ``H`` is written as a formal linear combination of
these edges and represented by a vector in ``\mathbb{Z}^{2d}``, where ``\sigma_i`` corresponds to the ``i``-th and
``\xi_i`` to the ``(i + d)``-th standard basis vector. This representation is not unique, since ``H`` is a subquotient
of ``\mathbb{Z}^{2d}``.

It is assumed that `o` has no nontrivial translations.

# Examples
```jldoctest
julia> o = origami(cperm([1,2,3,4,5],[6,7]), cperm([1,6,8],[2,7]))
Origami ((1,2,3,4,5)(6,7),(1,6,8)(2,7), 8)

julia> homology(o)
6-element Vector{Any}:
 GAP: [ 0, 0, 1, 1, 1, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0 ]
 GAP: [ 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, -2, 0 ]
 GAP: [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0 ]
 GAP: [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ]
 GAP: [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ]
 GAP: [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ]
```
"""
function homology(o::Origami)
  return GAP.gap_to_julia(GAP.Globals.HomologyOrigami(GapObj(o)))
end

@doc raw"""
    non_taut_part_of_homology(o::Origami, homology_basis)

Return a ``\mathbb{Z}``-basis of the non-tautological part ``H^{NT}`` of the homology of the origami `o`, with respect
to the basis `homology_basis` of the homology, for example the one returned by [`homology`](@ref).

``H^{NT}`` is the kernel of ``p_*\colon H \to H_1(E, \mathbb{Z}) \cong \mathbb{Z}^2`` for the covering ``p\colon X \to E``
of the torus. It has rank ``2g - 2``, where ``g`` is the genus of `o`, is invariant under the action of the affine group
of `o`, and is the orthogonal complement of the tautological part ``H^T`` of rank ``2`` with respect to the intersection
form. Over ``\mathbb{Q}`` the homology is the direct sum of ``H^T`` and ``H^{NT}``; over ``\mathbb{Z}``,
``H^T \oplus H^{NT}`` is in general only a subgroup of finite index in ``H``.

It is assumed that `o` has no nontrivial translations.

# Examples
```jldoctest
julia> o = origami(cperm([1,2,3,4,5],[6,7]), cperm([1,6,8],[2,7]))
Origami ((1,2,3,4,5)(6,7),(1,6,8)(2,7), 8)

julia> non_taut_part_of_homology(o, homology(o))
4-element Vector{Any}:
 GAP: [ 0, 0, 1, 1, 1, -1, 0, -2, 1, 0, 2, -1, 0, 0, -1, -1 ]
 GAP: [ 0, 0, 0, 0, 0, 0, 1, -1, 0, -1, 1, 0, 0, 0, -2, 2 ]
 GAP: [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 ]
 GAP: [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1 ]
```
"""
function non_taut_part_of_homology(o::Origami, homology_basis)
  return GAP.gap_to_julia(
    GAP.Globals.NonTautPartOfHomologyOrigami(GapObj(o), GapObj(homology_basis)))
end

@doc raw"""
    action_of_matrix_on_homology(o::Origami, a::Matrix)

Return the action of an affine homeomorphism of the origami `o` with derivative `a` on the homology of `o`, with respect
to the basis returned by [`homology`](@ref). The result is a matrix in ``{\rm Sp}_{2g}(\mathbb{Z})``, given as the list
of its rows.

It is assumed that `o` has no nontrivial translations.

# Examples
```jldoctest
julia> o = origami(cperm([1,2,3,4,5],[6,7]), cperm([1,6,8],[2,7]))
Origami ((1,2,3,4,5)(6,7),(1,6,8)(2,7), 8)

julia> action_of_matrix_on_homology(o, [1 10; 0 1])
6-element Vector{Any}:
 GAP: [ 1, -2, -2, 2, 0, 0 ]
 GAP: [ -5, -11, -2, 2, 5, 0 ]
 GAP: [ -5, -14, -3, 4, 5, 10 ]
 GAP: [ -10, -26, -6, 7, 10, 10 ]
 GAP: [ -10, -26, -6, 6, 11, 0 ]
 GAP: [ 0, 0, 0, 0, 0, 1 ]
```
"""
function action_of_matrix_on_homology(o::Origami, a::Matrix)
  return GAP.gap_to_julia(
    GAP.Globals.ActionOfMatrixOnHom(GapObj(o), GapObj(a)))
end

@doc raw"""
    action_of_matrix_on_non_taut(o::Origami, a::Matrix)

Return the action of an affine homeomorphism of the origami `o` with derivative `a` on the non-tautological part of the
homology of `o`. The bases are the ones returned by [`homology`](@ref) and [`non_taut_part_of_homology`](@ref). The
result is a matrix in ``{\rm Sp}_{2g-2}(\mathbb{Z})``, given as the list of its rows.

It is assumed that `o` has no nontrivial translations.

# Examples
```jldoctest
julia> o = origami(cperm([1,2,3,4,5],[6,7]), cperm([1,6,8],[2,7]))
Origami ((1,2,3,4,5)(6,7),(1,6,8)(2,7), 8)

julia> Oscar.action_of_matrix_on_non_taut(o, [1 10; 0 1])
4-element Vector{Any}:
 GAP: [ 5, 0, 2, 0 ]
 GAP: [ -1, -9, 2, 5 ]
 GAP: [ -8, 0, -3, 0 ]
 GAP: [ 2, -20, 6, 11 ]
```
"""
function action_of_matrix_on_non_taut(o::Origami, a::Matrix)
  return GAP.gap_to_julia(
    GAP.Globals.ActionOfMatrixOnNonTaut(GapObj(o), GapObj(a)))
end

@doc raw"""
    shadow_veech_group(o::Origami)

Return the shadow Veech group of the origami `o`, that is, the image of the action of the Veech group of `o` on the
non-tautological part of the homology of `o`. It is a matrix group contained in ``{\rm Sp}_{2g-2}(\mathbb{Z})``. The
bases are the ones returned by [`homology`](@ref) and [`non_taut_part_of_homology`](@ref).

It is assumed that `o` has no nontrivial translations.

# Examples
```jldoctest
julia> o = origami(cperm([1,2,3,4,5],[6,7]), cperm([1,6,8],[2,7]))
Origami ((1,2,3,4,5)(6,7),(1,6,8)(2,7), 8)

julia> Oscar.shadow_veech_group(o)
Matrix group of degree 4
  over integer ring
```
"""
function shadow_veech_group(o::Origami)
  shadowVeech = GAP.Globals.ShadowVeechGroup(GapObj(o))::GapObj
  generators = GAP.Globals.GeneratorsOfGroup(shadowVeech)::GapObj
  julia_matrices = [matrix(ZZ, gen) for gen in generators]
  return matrix_group(julia_matrices)
end

function homology_to_string(h)
  # TODO
end
