using PTpost, Test, Infiltrator
import Plasmatrace

function InternochronTest()
    Plasmatrace.PT(PTpost,logbook="test.log")
end

@testset "Internochron test" begin InternochronTest() end
