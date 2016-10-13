@testset "identities and expectional cases for $T" for T in (Float64, Float32,)

    @test isnan(Amal.exp(T(NaN)))
    @test Amal.exp(T(-Inf)) === T(0)
    @test Amal.exp(T(Inf)) === T(Inf)
    @test Amal.exp(T(0)) === T(1) # exact

    for x in T[10000, -10000]
        @test cmpdenorm(T, Amal.exp(x), Base.exp(BigFloat(x)))
    end

    @test isnan(Amal.exp2(T(NaN)))
    @test Amal.exp2(T(-Inf)) === T(0)
    @test Amal.exp2(T(Inf)) === T(Inf)
    @test Amal.exp2(T(0)) === T(1) # exact
    @test Amal.exp2(T(2)) === T(4)
    @test Amal.exp2(T(12)) === T(4096)

    for x in T[10000, -10000]
        @test cmpdenorm(T, Amal.exp2(x), Base.exp2(BigFloat(x)))
    end

    @test isnan(Amal.exp10(T(NaN)))
    @test Amal.exp10(T(-Inf)) === T(0)
    @test Amal.exp10(T(Inf)) === T(Inf)
    @test Amal.exp10(T(0)) === T(1) # exact
    @test Amal.exp10(T(1)) === T(10)
    @test Amal.exp10(T(3)) === T(1000)

    for x in T[10000, -10000]
        @test cmpdenorm(T, Amal.exp10(x), Base.exp10(BigFloat(x)))
    end

end
