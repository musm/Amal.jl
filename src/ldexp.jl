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
    ldexp(a::IEEEFloat, n::Int) -> Float

Computes `a × 2^n`.
"""
@inline function _ldexp{T<:IEEEFloat}(x::T, q::Integer)
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
