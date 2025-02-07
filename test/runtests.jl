using PileResponse
using Test

@testset "PileResponse.jl" begin
    @test PileResponse.helloworld() == 123
end
