if !(@isdefined rerun)
    using Revise, Pkg
    Pkg.activate("/home/pvermees/git/PTpost.jl")
    Pkg.instantiate()
    Pkg.precompile()
    cd("/home/pvermees/git/PTpost.jl/test")
end

rerun = true

include("runtests.jl")
