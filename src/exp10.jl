#  Method
#    1. Argument reduction: Reduce x to an r so that |r| <= 0.5*log10(2). Given x,
#       find r and integer k such that
#
#                x = k*log10(2) + r,  |r| <= 0.5*log10(2).
#
#    2. Approximate exp10(r) by a polynomial on the interval [-0.5*log10(2), 0.5*log10(2)]:
#
#           exp10(x) = 1.0 + polynomial(x),
#
#    3. Scale back: exp10(x) = 2^k * exp10(r)

@inline exp10_kernel{T<:LargeFloat}(x::T) =  @horner_oftype(x, 1.0,
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

@inline exp10_kernel{T<:SmallFloat}(x::T) = @horner_oftype(x, 1.0,
    2.302585124969482421875,
    2.650949001312255859375,
    2.0346698760986328125,
    1.17125606536865234375,
    0.5400512218475341796875,
    0.20749187469482421875,
    5.2789829671382904052734375e-2)

"""
    exp10(x)

Compute the base `10` exponential of `x`, in other words ``10^x``.
"""
function exp10{T<:IEEEFloat}(x::T)
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
        x > MAXEXP10(T) && return T(Inf)
        x < MINEXP10(T) && return T(0.0)
    end
    
    # reduce
    k = round(T(LOG210)*x)
    n = _trunc(k)
    r = muladd(k, -LOG102U(T), x)
    r = muladd(k, -LOG102L(T), r)

    # compute approximation
    u = exp10_kernel(r)
    return _ldexp(u, n)
end
