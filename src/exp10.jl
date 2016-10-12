"""
    exp10(x)

Compute the base ``10`` exponential of ``x``, in other words ``10^x``.
"""
function exp10 end

# Method: see exp.jl, same idea except we do not have the hi, lo argument split

@inline @oftype_float function _exp10{T}(r::T)
    z = r*r
    p = 0.86858896380650363333586483349790796637535095214844 + z *
    (0.38376418216572089070481865746842231601476669311523 + z *
    (-3.3911309900390991523000394636255805380642414093018e-2 + z *
    (4.2808249045911085303717236172360571799799799919128e-3 + z *
    (-5.6859012685279270403471141293039181618951261043549e-4 + z *
    (1.7440650302036995599082314090111367477220483124256e-4 + z *
    (-4.6916542309206350075401203980618447531014680862427e-3 + z *
    (0.118128643884600606495105523663369240239262580871582 + z *
    (-1.2272881009328624468679436176898889243602752685547))))))))
    return 1.0 + 2.0*r/(p - r)
end

# @inline @oftype_float function _exp10{T<:SmallFloatTypes}(r::T)
#     z = r*r
#     p = 0.868588924407958984375 + z * 
#     (0.383812963962554931640625 + z * 
#     (-5.1157988607883453369140625e-2 + z *
#     (2.616692066192626953125 + z * 
#     (-192.3223876953125 + z * 
#     (6795.14453125 + z *
#     (-92332.234375))))))
#     return 1.0 + 2.0*r/(p - r)
# end

@oftype_float function exp10{T}(x::T)
    x > MAXEXP10(T) && return Inf
    x < MINEXP10(T) && return 0.0
 
    # reduce
    k = round(T(LOG210)*x) 
    n = _trunc(k)
    r = muladd(k, -LOG102U(T), x)
    r = muladd(k, -LOG102L(T), r)

    # compute approximation
    x = _exp10(r)
    return _ldexp(x,n)
end