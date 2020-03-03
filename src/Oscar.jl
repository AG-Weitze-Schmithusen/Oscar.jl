module Oscar

#=
  We currently only import packages which:
    * are registered
    * have BinaryBuilder or build in a few minutes
    * support Linux and OSX
=#

import AbstractAlgebra
import Nemo
import Hecke
import Singular
import Polymake
using Markdown

import AbstractAlgebra: @show_name, @show_special

function __init__()
  println(" -----    -----    -----      -      -----   ")
  println("|     |  |     |  |     |    | |    |     |  ")
  println("|     |  |        |         |   |   |     |  ")
  println("|     |   -----   |        |     |  |-----   ")
  println("|     |        |  |        |-----|  |   |    ")
  println("|     |  |     |  |     |  |     |  |    |   ")
  println(" -----    -----    -----   -     -  -     -  ")
  println()
  println("...combining (and extending) Gap, Hecke, Nemo, Polymake and Singular")
end

#function load()
#  s = joinpath(dirname(pathof(Oscar)), "Modules/FreeModules-graded.jl")
#  Base.include(Oscar, s)
#  Revise.track(Oscar, s)
#end

include("OscarTypes.jl")

include("Rings/integer.jl")
include("Rings/rational.jl")
include("Rings/Hecke.jl")
include("Rings/mpoly.jl")
include("Rings/mpoly-graded.jl")
include("Modules/FreeModules-graded.jl")
include("Polymake/Ineq.jl")
include("Polymake/NmbThy.jl")

end # module
