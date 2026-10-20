# Build only the Origami part of the OSCAR documentation, which is much faster
# than `Oscar.build_doc()`, and open it in the browser. Doctests are not run.
#
# Usage in the Julia REPL (`Oscar.build_origami_doc` loads this file on demand):
#   using Revise, Oscar
#   Oscar.build_origami_doc()
#
# Load Revise before Oscar, so that docstrings edited afterwards are picked up
# by further calls without restarting Julia.
#
# Usage from the shell:
#   julia --project=. experimental/Origami/docs/build_origami_doc.jl [outdir]
module OrigamiDocBuilder

using Oscar

# Set up the documentation environment only once per session; this makes
# Documenter loadable without switching the active project.
isdefined(Oscar, :docsproject) || Oscar.doc_init()
using Documenter, DocumenterCitations

"""
    build_origami_doc(outdir::String = joinpath(tempdir(), "origami_doc"); open_browser::Bool = true)

Build the Origami documentation pages into `outdir` and, if `open_browser`
is `true`, open them in the browser. `outdir` is deleted and recreated on
every call. Return the path of the generated `origamis.html`.
"""
function build_origami_doc(outdir::String = joinpath(tempdir(), "origami_doc"); open_browser::Bool = true)
  outdir = abspath(outdir)
  oscardir = Base.pkgdir(Oscar)
  srcdir = joinpath(outdir, "src")
  rm(outdir; recursive=true, force=true)
  mkpath(srcdir)

  for f in ("introduction.md", "origamis.md")
    cp(joinpath(oscardir, "experimental", "Origami", "docs", "src", f), joinpath(srcdir, f))
  end
  write(joinpath(srcdir, "references.md"), "# Bibliography\n\n```@bibliography\n```\n")

  DocMeta.setdocmeta!(Oscar, :DocTestSetup, Oscar.doctestsetup(); recursive=true, warn=false)
  bib = CitationBibliography(joinpath(oscardir, "docs", "oscar_references.bib"))

  Oscar.with_unicode(false) do
    makedocs(;
      root=outdir,
      source="src",
      build="build",
      format=Documenter.HTML(; prettyurls=false),
      sitename="Origami (OSCAR)",
      modules=[Oscar],
      doctest=false,
      meta=Oscar.docmeta(),
      warnonly=true,
      checkdocs=:none,
      remotes=nothing,
      pages=["Origami" => ["introduction.md", "origamis.md"], "references.md"],
      plugins=[bib],
    )
  end

  page = joinpath(outdir, "build", "origamis.html")
  println("\nDocumentation built: ", page)
  open_browser && _open_in_browser(page)
  return page
end

function _open_in_browser(file::String)
  if Sys.isapple()
    run(`open $file`; wait=false)
  elseif Sys.iswindows()
    run(`$(get(ENV, "COMSPEC", "cmd.exe")) /c start $file`; wait=false)
  elseif !isnothing(Sys.which("wslview")) # Linux under WSL
    run(`wslview $file`; wait=false)
  elseif !isnothing(Sys.which("xdg-open"))
    run(`xdg-open $file`; wait=false)
  else
    @warn "Could not open a browser, please open $file manually."
  end
end

end # module OrigamiDocBuilder

# Build directly when run as a script from the shell.
if abspath(PROGRAM_FILE) == @__FILE__
  OrigamiDocBuilder.build_origami_doc(isempty(ARGS) ? joinpath(tempdir(), "origami_doc") : ARGS[1])
end
