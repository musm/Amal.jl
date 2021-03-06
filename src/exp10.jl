#  Method
#  1. Argument reduction: Reduce x to an r so that |r| <= 0.5*log10(2). Given x,
#     find r and integer k such that
#
#                x = k*log10(2) + r,  |r| <= 0.5*log10(2).
#
# 2. Approximate exp10(r) by a polynomial on the interval [-0.5*log10(2), 0.5*log10(2)]:
#
#           exp10(x) = 1.0 + polynomial(x),
#
#    sup norm relative error within the interval of the polynomial approximations:
#    Float64 : [2.7245504724394698952e-18; 2.7245529895753476720e-18]
#    Float32 : [9.6026471477842205871e-10; 9.6026560194009888672e-10]
#
# 3. Scale back: exp10(x) = 2^k * exp10(r)

@inline exp10_kernel(x::Float64) =
    @horner(x, 1.0,
    2.30258509299404590109361379290930926799774169921875,
    2.6509490552391992146397114993305876851081848144531,
    2.03467859229323178027470930828712880611419677734375,
    1.17125514891212478829629617393948137760162353515625,
    0.53938292928868392106522833273629657924175262451172,
    0.20699584873167015119932443667494226247072219848633,
    6.8089348259156870502017966373387025669217109680176e-2,
    1.9597690535095281527677713029333972372114658355713e-2,
    5.015553121397981796436571499953060992993414402008e-3,
    1.15474960721768829356725927226534622604958713054657e-3,
    1.55440426715227567738830671828509366605430841445923e-4,
    3.8731032432074128681303432086835414338565897196531e-5,
    2.3804466459036747669197886523306806338950991630554e-3,
    9.3881392238209649520573607528461934634833596646786e-5,
    -2.64330486232183387018679354696359951049089431762695e-2)

@inline exp10_kernel(x::Float32) =
    @horner(x, 1.0f0,
    2.302585124969482421875f0,
    2.650949001312255859375f0,
    2.0346698760986328125f0,
    1.17125606536865234375f0,
    0.5400512218475341796875f0,
    0.20749187469482421875f0,
    5.2789829671382904052734375f-2)

@eval exp10_small_thres(::Type{Float64}) = $(2.0^-29)
@eval exp10_small_thres(::Type{Float32}) = $(2.0f0^-14)

"""
    exp10(x)

Compute the base `10` exponential of `x`, in other words ``10^x``.

```jldoctest
julia> exp10(1.0)
10.0
```
"""
function exp10(x::T) where {T<:Union{Float32,Float64}}
    xa = reinterpret(Unsigned, x) & ~sign_mask(T)
    xsb = signbit(x)

    # filter out non-finite arguments
    if xa > reinterpret(Unsigned, EXP10_MAX(T))
        if xa >= exponent_mask(T)
            xa & significand_mask(T) != 0 && return T(NaN)
            return xsb ? T(0.0) : T(Inf) # exp10(+-Inf)
        end
        x > EXP10_MAX(T) && return T(Inf)
        x < EXP10_MIN(T) && return T(0.0)
    end
    # compute approximation
    if xa > reinterpret(Unsigned, T(0.5)*T(LOG10_2)) # |r| > 0.5 log10(2).
        # argument reduction
        if xa < reinterpret(Unsigned, T(1.5)*T(LOG10_2)) # |r| <= 1.5 log10(2)
            if xsb
                k = -1
                r = LOG10_2U(T) + x
                r = LOG10_2L(T) + r
            else
                k = 1
                r = x - LOG10_2U(T)
                r = r - LOG10_2L(T)
            end
        else
            n = round(T(LOG2_10)*x)
            k = unsafe_trunc(Int,n)
            r = muladd(n, -LOG10_2U(T), x)
            r = muladd(n, -LOG10_2L(T), r)
        end
        # compute approximation on reduced argument
        y = exp10_kernel(r)
        # scale back
        if k > -significand_bits(T)
            # multiply by 2.0 first to prevent overflow, extending the range
            k == exponent_max(T) && return y * T(2.0) * T(2.0)^(exponent_max(T) - 1)
            twopk = reinterpret(T, rem(exponent_bias(T) + k, fpinttype(T)) << significand_bits(T))
            return y*twopk
        else
            # add significand_bits(T) + 1 to lift the range outside the subnormals
            twopk = reinterpret(T, rem(exponent_bias(T) + significand_bits(T) + 1 + k, fpinttype(T)) << significand_bits(T))
            return y * twopk * T(2.0)^(-significand_bits(T) - 1)
        end
    elseif xa < reinterpret(Unsigned, exp10_small_thres(T))  # |x| < exp10_small_thres
        # Taylor approximation for small values: exp10(x) ≈ 1.0 + log(10)*x
        return muladd(x, T(LN10), T(1.0))
    else
        # primary range with k = 0, so compute approximation directly
        return exp10_kernel(x)
    end
end
