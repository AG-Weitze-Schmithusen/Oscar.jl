@doc raw"""
    build_origami_doc(outdir::String = joinpath(tempdir(), "origami_doc"); open_browser::Bool = true)

Build only the Origami part of the OSCAR documentation into `outdir` and, if
`open_browser` is `true`, open it in the browser. This is much faster than
[`build_doc`](@ref), but does not run doctests. `outdir` is deleted and
recreated on every call. Return the path of the generated `origamis.html`.

Load Revise before Oscar, so that docstrings edited afterwards are picked up
by further calls without restarting Julia.
"""
function build_origami_doc(outdir::String = joinpath(tempdir(), "origami_doc"); open_browser::Bool = true)
  # Load the builder only on demand, so that Oscar does not depend on Documenter.
  if !isdefined(Main, :OrigamiDocBuilder)
    Base.include(Main, joinpath(@__DIR__, "..", "docs", "build_origami_doc.jl"))
  end
  return Base.invokelatest(() -> Main.OrigamiDocBuilder.build_origami_doc(outdir; open_browser))
end
