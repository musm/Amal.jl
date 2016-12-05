"""
    ldexp(x::AbstractFloat, n::Integer)

Computes ``x \\times 2^n``.
"""
function ldexp{T<:AbstractFloat}(x::T, e::Integer)
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xs >= exponent_mask(T) && return x # NaN or Inf
    k = Int(xs >> significand_bits(T))
    if xs <= (~exponent_mask(T) & ~sign_mask(T)) # x is subnormal
        xs == 0 && return x # +-0
        m = leading_zeros(xs) - exponent_bits(T)
        ys = xs << unsigned(m)
        xu = ys | (xu & sign_mask(T))
        k = 1 - m
        # underflow, otherwise may have integer underflow in the following n + k
        e < -50000 && return flipsign(T(0.0), x)
    end
    # for cases larger than Int make sure we properly overlfow/underflow (optimized away)
    if e > typemax(Int)
        return flipsign(T(Inf), x)
    elseif e < typemin(Int)
        return flipsign(T(0.0), x)
    end
    n = e % Int
    k += n
    # overflow, if k is larger than maximum posible exponent
    if k >= exponent_raw_max(T)
        return flipsign(T(Inf), x)
    end
    if k > 0 # normal case
        xu = (xu & ~exponent_mask(T)) | (k % typeof(xu) << significand_bits(T))
        return reinterpret(T, xu)
    else # subnormal case
        if k <= -significand_bits(T) # underflow
            # overflow, for the case of integer overflow in n + k
            e > 50000 && return flipsign(T(Inf), x)
            return flipsign(T(0.0), x)
        end
        k += significand_bits(T)
        z = T(2.0)^-significand_bits(T)
        xu = (xu & ~exponent_mask(T)) | (k % typeof(xu) << significand_bits(T))
        return z*reinterpret(T, xu)
    end
end


# private

@inline split_exponent(::Type{Float64}, q::Int) = _split_exponent(q, UInt(9), UInt(31), UInt(2))
@inline split_exponent(::Type{Float32}, q::Int) = _split_exponent(q, UInt(6), UInt(31), UInt(2))
@inline split_exponent(::Type{Float16}, q::Int) = _split_exponent(q, UInt(3), UInt(31), UInt(2))
@inline function _split_exponent(q, n, v, offset)
    m = q >> v
    m = (((m + q) >> n) - m) << (n-offset)
    q = q - (m << offset)
    return m, q
end

intasfloat{T<:AbstractFloat}(::Type{T}, m::Integer) = reinterpret(T, (inttype(T)(m)) << significand_bits(T))

@inline function _ldexp{T<:AbstractFloat}(x::T, q::Integer)
    bias = exponent_bias(T)
    emax = exponent_raw_max(T)
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
