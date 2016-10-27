inttype(::Type{Float64}) = Int64
inttype(::Type{Float32}) = Int32
asint{T}(::Type{T}, x) = inttype(T)(x)

function NaNs{T}(::Type{T}, i) # i starts counting from the most significant bit
    i += 1
    reinterpret(T,(asint(T, Base.exponent_mask(T)) | asint(T,1) << (Base.significand_bits(T) - i))) # most sig bit
end

@testset "identities and expectional cases for $T" for T in (Float64, Float32)

    @test isnan(Amal.frexp(T(NaN))[1]) && Amal.frexp(T(NaN))[2] == 0
    @test isnan(Amal.frexp(NaNs(T,1))[1]) && Amal.frexp(NaNs(T,1))[2] == 0
    @test isnan(Amal.frexp(NaNs(T,5))[1]) && Amal.frexp(NaNs(T,5))[2] == 0
    @test Amal.frexp(T(0.0)) === (T(0.0), 0)
    @test Amal.frexp(T(-0.0)) === (T(-0.0), 0)
    @test Amal.frexp(T(Inf)) == (T(Inf), 0)
    @test Amal.frexp(T(-Inf)) == (T(-Inf), 0)

    @test isnan(Amal.exp(T(NaN)))
    @test isnan(Amal.exp(NaNs(T,1)))
    @test isnan(Amal.exp(NaNs(T,5)))
    @test Amal.exp(T(-Inf))      == T(0.0)
    @test Amal.exp(T(Inf))       == T(Inf)
    @test Amal.exp(T(0.0))       == T(1.0) # exact
    @test Amal.exp(T(5000.0))    == T(Inf)
    @test Amal.exp(T(-5000.0))   == T(0.0)

    @test isnan(Amal.exp2(T(NaN)))
    @test isnan(Amal.exp2(NaNs(T,1)))
    @test isnan(Amal.exp2(NaNs(T,5)))
    @test Amal.exp2(T(-Inf))     == T(0.0)
    @test Amal.exp2(T(Inf))      == T(Inf)
    @test Amal.exp2(T(0.0))      == T(1.0) # exact
    @test Amal.exp2(T(2.0))      == T(4.0)
    @test Amal.exp2(T(12.0))     == T(4096.0)
    @test Amal.exp2(T(5000.0))   == T(Inf)
    @test Amal.exp2(T(-5000.0))  == T(0.0)

    @test isnan(Amal.exp10(T(NaN)))
    @test isnan(Amal.exp10(NaNs(T,1)))
    @test isnan(Amal.exp10(NaNs(T,5)))
    @test Amal.exp10(T(-Inf))    == T(0.0)
    @test Amal.exp10(T(Inf))     == T(Inf)
    @test Amal.exp10(T(0.0))     == T(1.0) # exact
    @test Amal.exp10(T(1.0))     == T(10.0)
    @test Amal.exp10(T(3.0))     == T(1000.0)
    @test Amal.exp10(T(5000.0))  == T(Inf)
    @test Amal.exp10(T(-5000.0)) == T(0.0)

    @test isnan(Amal.log(T(NaN)))
    @test isnan(Amal.log(NaNs(T,1)))
    @test isnan(Amal.log(NaNs(T,5)))
    @test Amal.log(T(Inf))       == T(Inf)
    @test Amal.log(T(0.0))       == T(-Inf)
    @test Amal.log(T(-0.0))      == T(-Inf)
    @test Amal.log(T(1.0))       == T(0.0) # exact
    @test_throws DomainError Amal.log(T(-1.0))

    # @test isnan(Amal.log2(T(NaN)))
    # @test Amal.log2(T(Inf))     == T(Inf)
    # @test Amal.log2(T(0.0))     == T(-Inf)
    # @test Amal.log2(T(-0.0))    == T(-Inf)
    # @test Amal.log2(T(1.0))     == T(0.0) # exact
    # @test_throws DomainError Amal.log2(T(-1.0))

    # @test isnan(Amal.cbrt(T(NaN)))
    # @test Amal.cbrt(T(-Inf)) == T(-Inf)
    # @test Amal.cbrt(T(Inf))  == T(Inf)
    # @test Amal.cbrt(T(0))    == = T(0)
    # @test Amal.cbrt(T(-0))   == = T(-0)
end
