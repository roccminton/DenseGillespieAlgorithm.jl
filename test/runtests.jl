using DenseGillespieAlgorithm
using Test

@testset "sumsumdict function" begin
    @test DenseGillespieAlgorithm.sumsumdict(Dict(
        "a" => [1,2,3],
        "b" => [4,5,6],
        "c" => [7,8,9]
    )) == 45
end

@testset "historylength function" begin
    #test length
    l = 100
    #default parameter
    par = (a=1,)

    #test vector
    v = zeros(l)
    @test DenseGillespieAlgorithm.historylength(v,par) == l

    #test matrix
    m = zeros(l,3)
    @test DenseGillespieAlgorithm.historylength(m,par) == l

    #test dict
    d = Dict("a" => zeros(l))
    @test DenseGillespieAlgorithm.historylength(d,par) == l

    #test par
    p = (historylength = l,)
    @test DenseGillespieAlgorithm.historylength("something",p) == l
end

@testset "dropzerso! function" begin
    #Test Number 1
    d = Dict("a" => 1.0, "b" => 0.0)
    expected = Dict("a" => 1.0)
    DenseGillespieAlgorithm.dropzeros!(d)
    @test d == expected

    #Test Number 2
    d = Dict("a" => 1.0, "b" => 0.0, c => 0.0)
    expected = Dict("a" => 1.0)
    DenseGillespieAlgorithm.dropzeros!(d)
    @test d == expected

    #Test Number 3
    d = Dict("a" => 1.0, "b" => 2.0, c => 3.0)
    expected = Dict("a" => 1.0, "b" => 2.0, c => 3.0)
    DenseGillespieAlgorithm.dropzeros!(d)
    @test d == expected

    #Test Vector 1
    d = Dict("a" => [1,2,3],"b" => [0,2,3])
    expected = Dict("a" => [1,2,3])
    DenseGillespieAlgorithm.dropzeros!(d)
    @test d == expected

    #Test Vector 2
    d = Dict("a" => [1,2,3],"b" => [0,2,3],"c" => [0,4,5])
    expected = Dict("a" => [1,2,3])
    DenseGillespieAlgorithm.dropzeros!(d)
    @test d == expected

    #Test Vector 3
    d = Dict("a" => [1,2,3],"b" => [1,2,3])
    expected = Dict("a" => [1,2,3],"b" => [1,2,3])
    DenseGillespieAlgorithm.dropzeros!(d)
    @test d == expected

    #Test Vector 3
    d = Dict("a" => [1,2,3],"b" => [1,0,3])
    expected = Dict("a" => [1,2,3],"b" => [1,0,3])
    DenseGillespieAlgorithm.dropzeros!(d)
    @test d == expected

    #Test other types 1
    d = [1,2,3]
    expected = [1,2,3]
    DenseGillespieAlgorithm.dropzeros!(d)
    @test d == expected

    #Test other types 2
    d = [0,2,3]
    expected = [0,2,3]
    DenseGillespieAlgorithm.dropzeros!(d)
    @test d == expected
end
