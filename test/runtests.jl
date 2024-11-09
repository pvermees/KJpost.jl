using PTpost, Test, Infiltrator, Plots
import Plasmatrace

function InternochronTest()
    myrun = Plasmatrace.load("Lu-Hf",instrument="Agilent")
    method = "Lu-Hf"
    channels = Dict("d"=>"Hf178 -> 260",
                    "D"=>"Hf176 -> 258",
                    "P"=>"Lu175 -> 175")
    standards = Dict("Hogsbo" => "hogsbo")
    glass = Dict("NIST612" => "NIST612p")
    blk, fit = Plasmatrace.process!(myrun,method,channels,standards,glass,
                                    nblank=2,ndrift=2,ndown=2)
    P, D, d = Plasmatrace.atomic(myrun[1],channels,blk,fit)
    x0, sx0, y0, sy0, rho = internochron(P,D,d)
    p = PTpost.plot(x0,y0,P,D,d)
    @test display(p) != NaN
end

@testset "Internochron test" begin InternochronTest() end
