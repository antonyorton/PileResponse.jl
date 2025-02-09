using PileResponse
using Test

@testset "PileResponse.jl" begin
    @test PileResponse.helloworld() == 123
    @test square(5) == 25
    @test square(4.0) == 16.0
end
