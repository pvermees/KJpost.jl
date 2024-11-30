if !(@isdefined rerun)
    using Revise, Pkg
    Pkg.activate("/home/pvermees/git/KJpost.jl")
    Pkg.instantiate()
    Pkg.precompile()
    cd("/home/pvermees/git/KJpost.jl/test")
end

rerun = true

include("runtests.jl")
