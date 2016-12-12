#  Method
#    1. Argument reduction: Reduce x to an r so that |r| <= 0.5. Given x,
#       find r and integer k such that
#
#                x = k + r,  |r| <= 0.5.
#
#    2. Approximate exp2(r) by a polynomial on the interval [-0.5, 0.5]:
#
#           exp2(x) = 1.0 + polynomial(x),
#
#    3. Scale back: exp2(x) = 2^k * exp2(r)

@inline exp2_kernel{T<:LargeFloat}(x::T) = @horner_oftype(x, 1.0,
    0.69314718055994528622676398299518041312694549560547,
    0.240226506959100721827482516346208285540342330932617,
    5.5504108664823380292485666132051846943795680999756e-2,
    9.618129107628051871481389412110729608684778213501e-3,
    1.3333558145962711230514408100589207606390118598938e-3,
    1.54035303940471083325794432461464111838722601532936e-4,
    1.52527343581534264687878457711356361414800630882382e-5,
    1.3215486321169555424209860611250988426945696119219e-6,
    1.01777526742633583156634826966807638726209006563295e-7,
    7.0550611328426153816964051335245550200525599393586e-9,
    4.5441884834923169985517595445740929999134394279281e-10,
    2.55208683164155542627190035541293782264671285986424e-11,
    -1.00458805873217110125204978507248889784547740688936e-11)

@inline exp2_kernel{T<:SmallFloat}(x::T) = @horner_oftype(x, 1.0,
    0.693147182464599609375,
    0.2402265071868896484375,
    5.5504061281681060791015625e-2,
    9.6180774271488189697265625e-3,
    1.333682797849178314208984375e-3,
    1.54563575051724910736083984375e-4,
    1.4599625501432456076145172119140625e-5)

"""
    exp2(x)

Compute the base `2` exponential of `x`, in other words ``2^x``.
"""
function exp2{T<:IEEEFloat}(x::T)
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xsb = xu & sign_mask(T)

    # filter out non-finite arguments
    if xs > reinterpret(Unsigned, MAXEXP(T))
        if xs >= exponent_mask(T)
            if xs & significand_mask(T) != 0
                return T(NaN) 
            end
            return xsb == 0 ? T(Inf) : T(0.0) # exp(+-Inf)
        end
        x > MAXEXP2(T) && return T(Inf)
        x < MINEXP2(T) && return T(0.0)
    end
    
    # reduce
    k = round(x)
    n = _trunc(k)
    r = x - k

    # compute approximation
    u = exp2_kernel(r)
    return _ldexp(u, n)
end
