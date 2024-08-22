module HFDisorderHubbard

export CubicLattice
export HubbardPara
export SCFdata
export init, step!

include("lattices.jl")
include("Hubbard.jl")
include("SCF.jl")

end
