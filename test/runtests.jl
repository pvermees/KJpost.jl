using Plasmatrace, PTpost, Test, Infiltrator

function PTguiTest()
    PT(PTpost;logbook="test.log")
end

function InternochronTest()
    myrun = load("Lu-Hf",instrument="Agilent")
    method = "Lu-Hf"
    channels = Dict("d"=>"Hf178 -> 260",
                    "D"=>"Hf176 -> 258",
                    "P"=>"Lu175 -> 175")
    standards = Dict("Hogsbo_gt" => "hogsbo")
    glass = Dict("NIST612" => "NIST612p")
    blk, fit = process!(myrun,method,channels,standards,glass,
                                    nblank=2,ndrift=2,ndown=2)
    P, D, d = atomic(myrun[1],channels,blk,fit)
    x0, y0, E = internochron(P,D,d)
    p = PTpost.plot(x0,y0,E,P,D,d;markercolor=:white)
    @test display(p) != NaN
end

@testset "PT test" begin PTguiTest() end
@testset "Internochron test" begin InternochronTest() end
