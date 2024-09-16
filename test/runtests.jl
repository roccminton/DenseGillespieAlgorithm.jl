using DenseGillespieAlgorithm
using Test

@testset "DenseGillespieAlgorithm.jl" begin
    @test DenseGillespieAlgorithm.sumsumdict(Dict(
        "a" => [1,2,3],
        "b" => [4,5,6],
        "c" => [7,8,9]
    )) == 45
end
