using Plasmatrace, PTpost, Test, Infiltrator

function PTguiTest()
    PT(PTpost;logbook="test.log")
end

function LuHfHelper()
    myrun = load("Lu-Hf",instrument="Agilent")
    method = "Lu-Hf"
    channels = Dict("d"=>"Hf178 -> 260",
                    "D"=>"Hf176 -> 258",
                    "P"=>"Lu175 -> 175")
    standards = Dict("Hogsbo_gt" => "hogsbo")
    glass = Dict("NIST612" => "NIST612p")
    blk, fit = process!(myrun,method,channels,standards,glass,
                        nblank=2,ndrift=2,ndown=2)
    return myrun, method, channels, blk, fit
end

function InternochronTest()
    myrun, method, channels, blk, fit = LuHfHelper()
    P, D, d = atomic(myrun[1],channels,blk,fit)
    x0, y0, E = internochron(P,D,d)
    p = PTpost.plot(x0,y0,E,P,D,d,method;markercolor=:white)
    @test display(p) != NaN
end

function UPbHelper()
    method = "U-Pb"
    myrun = load("carbonate",instrument="Agilent")
    standards = Dict("WC1_cc"=>"WC1")
    glass = Dict("NIST612"=>"NIST612")
    channels = Dict("d"=>"Pb207","D"=>"Pb206","P"=>"U238")
    blk, fit = process!(myrun,method,channels,standards,glass,
                        nblank=2,ndrift=1,ndown=1)
    return myrun, method, channels, blk, fit
end

function UPbTest()
    myrun, method, channels, blk, fit = UPbHelper()
    P, D, d = atomic(myrun[7],channels,blk,fit)
    x0, y0, E = internochron(P,D,d)
    p = PTpost.plot(x0,y0,E,P,D,d,method;
                    markercolor=:white)
    @test display(p) != NaN
end

function CSVtest()
    myrun, method, channels, blk, fit = UPbHelper()
    PTpost.internochron(myrun,method,channels,blk,fit)
end

@testset "PT test" begin PTguiTest() end
#@testset "Internochron test" begin InternochronTest() end
#@testset "UPb test" begin UPbTest() end
#@testset "CSV test" begin CSVtest() end
