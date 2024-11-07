using PTpost, Test, Infiltrator
import Plasmatrace

function ExtensionTest()
    Plasmatrace.PT(PTpost)#,logbook="test.log")
end

@testset "Extension test" begin ExtensionTest() end
