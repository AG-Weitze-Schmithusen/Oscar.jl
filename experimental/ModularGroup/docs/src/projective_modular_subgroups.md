# Projective modular subgroups

## Construction of projective modular subgroups

```@docs
projective_modular_subgroup_via_right_action
projective_modular_subgroup_via_left_action
s_right_perm(G::ProjectiveModularGroup)
t_right_perm(G::ProjectiveModularGroup)
r_right_perm(G::ProjectiveModularGroup)
j_right_perm(G::ProjectiveModularGroup)
```

## Generators and membership

```@docs
word_gens(G::ProjectiveModularGroup)
projective_s_t_presentation
is_word_element_of(w::FPGroupElem, G::ProjectiveModularGroup)
coset_right_action_of(A::ZZMatrix, G::ProjectiveModularGroup)
coset_left_action_of(A::ZZMatrix, G::ProjectiveModularGroup)
defines_projective_coset_action
```
