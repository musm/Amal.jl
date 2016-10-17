@inline split_exponent(::Type{Float64}, q::Int) = _split_exponent(q, UInt(9), UInt(31), UInt(2))
@inline split_exponent(::Type{Float32}, q::Int) = _split_exponent(q, UInt(6), UInt(31), UInt(2))
@inline split_exponent(::Type{Float16}, q::Int) = _split_exponent(q, UInt(3), UInt(31), UInt(2))
@inline function _split_exponent(q, n, v, offset)
    m = q >> v
    m = (((m + q) >> n) - m) << (n-offset)
    q = q - (m << offset)
    return m, q
end

"""
    -ldexp(a::IntrinsicFloats, n::Int) -> Float

Computes `a × 2^n`.
"""
@inline function _ldexp{T<:IntrinsicFloats}(x::T, q::Integer)
    bias = exponent_bias(T)
    emax = exponent_max(T)
    m, q = split_exponent(T,q)
    m += bias
    m = ifelse(m < 0, 0, m)
    m = ifelse(m > emax, emax, m)
    q += bias
    u = intasfloat(T,m)
    x = x*u*u*u*u
    u = intasfloat(T,q)
    return x*u
end


# The following define threshold values for `ilog2k`
threshold_exponent(::Type{Float64}) = 300
threshold_exponent(::Type{Float32}) = 64

"""
    ilog2(x::Float) -> Int
Returns the integral part of the logarithm of `|x|`, using 2 as base for the logarithm; in other
words this returns the binary exponent of `x` so that
    x = significand × 2^exponent
where `significand ∈ [0.5, 1)`.
"""
@inline function _ilog2{T<:IntrinsicFloats}(d::T)
    bias = exponent_bias(T)
    emax = exponent_max(T)
    thres = threshold_exponent(T)

    m = d < T(2.0)^-thres
    d = ifelse(m, d*T(2)^thres, d)
    q = floatasint(d) & emax
    q = ifelse(m, q - (thres + bias - 1), q - (bias - 1)) # subtract 1 since we want 2^q
end

using Base: exponent_half, exponent_mask, significand_bits, leading_zeros, exponent_bits, sign_mask
@inline function _frexp{T<:IntrinsicFloats}(x::T)
    xu = reinterpret(Unsigned,x)
    xe = xu & exponent_mask(T)
    k = Int(xe >> significand_bits(T))
    if xe == 0 # x is subnormal
        x == 0 && return x, 0
        xs = xu & sign_mask(T)
        xu $= xs
        m = leading_zeros(xu)-exponent_bits(T)
        xu <<= m
        xu $= xs
        k = 1-m
    elseif xe == exponent_mask(T) # NaN or Inf
        return x,0
    end
    k -= (exponent_bias(T)-1)
    xu = (xu & ~exponent_mask(T)) | exponent_half(T)
    reinterpret(T,xu), k
end