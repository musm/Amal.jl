@testset "identities and expectional cases for $T" for T in (Float64, Float32)

    @test isnan(Amal.exp(T(NaN)))
    @test Amal.exp(T(-Inf))    == T(0)
    @test Amal.exp(T(Inf))     == T(Inf)
    @test Amal.exp(T(0))       == T(1) # exact
    @test Amal.exp(T(5000))    == T(Inf)
    @test Amal.exp(T(-5000))   == T(0)

    @test isnan(Amal.exp2(T(NaN)))
    @test Amal.exp2(T(-Inf))   == T(0)
    @test Amal.exp2(T(Inf))    == T(Inf)
    @test Amal.exp2(T(0))      == T(1) # exact
    @test Amal.exp2(T(2))      == T(4)
    @test Amal.exp2(T(12))     == T(4096)
    @test Amal.exp2(T(5000))   == T(Inf)
    @test Amal.exp2(T(-5000))  == T(0)

    @test isnan(Amal.exp10(T(NaN)))
    @test Amal.exp10(T(-Inf))  == T(0)
    @test Amal.exp10(T(Inf))   == T(Inf)
    @test Amal.exp10(T(0))     == T(1) # exact
    @test Amal.exp10(T(1))     == T(10)
    @test Amal.exp10(T(3))     == T(1000)
    @test Amal.exp10(T(5000))  == T(Inf)
    @test Amal.exp10(T(-5000)) == T(0)

    @test isnan(Amal.log(T(NaN)))
    @test Amal.log(T(Inf))     == T(Inf)
    @test Amal.log(T(0))       == T(-Inf)
    @test Amal.log(T(-0))      == T(-Inf)
    @test Amal.log(T(1))       == T(0) # exact
    @test_throws DomainError Amal.log(T(-1))

    # @test isnan(Amal.cbrt(T(NaN)))
    # @test Amal.cbrt(T(-Inf)) == T(-Inf)
    # @test Amal.cbrt(T(Inf))  == T(Inf)
    # @test Amal.cbrt(T(0))    == = T(0)
    # @test Amal.cbrt(T(-0))   == = T(-0)
end
