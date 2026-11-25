# Modular subgroups

## Construction of modular subgroups

```@docs
modular_subgroup_via_right_action
modular_subgroup_via_left_action
s_right_perm(G::ModularGroup)
t_right_perm(G::ModularGroup)
r_right_perm(G::ModularGroup)
j_right_perm(G::ModularGroup)
```

## Generators and membership

```@docs
word_gens(G::ModularGroup)
s_t_decomposition
s_t_presentation
is_word_element_of(w::FPGroupElem, G::ModularGroup)
right_coset_representatives
left_coset_representatives
coset_right_action_of(A::ZZMatrix, G::ModularGroup)
coset_left_action_of(A::ZZMatrix, G::ModularGroup)
defines_coset_action_s_t
```
