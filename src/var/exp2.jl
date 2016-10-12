"""
    exp2(x)

Compute the base ``2`` exponential of ``x``, in other words ``2^x``.
"""
function exp2 end

# Method: see exp.jl, same idea except we do not have the hi, lo argument split in this case

@inline @oftype_float function _exp2{T}(r::T)
    z = r*r
    p = 2.88539008177792677400930188014172017574310302734375 + z *
    (0.115524530093332036817521668581321137025952339172363 + z *
    (-9.2506847818508888391803024475734673615079373121262e-4 + z *
    (1.05822000247694508344355884821297308917564805597067e-5 + z *
    (-1.2725311062591074468333673941344841296086087822914e-7 + z *
    (2.6544023099239800410890694517792101625452971802588e-9 + z *
    (-4.8117604480494130690177356432513100514825055142865e-9 + z *
    (1.0959102409244011630035248725575230954731864585483e-8 + z *
    (-1.03171025586816884494039390542108325377057553851046e-8))))))))
    return 1.0 + 2.0*r/(p - r)
end

@inline @oftype_float function _exp2{T<:SmallFloatTypes}(r::T)
    z = r*r
    p = 2.88539028167724609375 + z * 
    (0.11550851166248321533203125 + z *
    (-5.65846334211528301239013671875e-4 + z *
    (-3.276503644883632659912109375e-3 + z *
    (1.315677352249622344970703125e-2 + z *
    (-1.91472284495830535888671875e-2)))))
    return 1.0 + 2.0*r/(p - r)
end

@oftype_float function exp2{T}(x::T)
    x > MAXEXP2(T) && return Inf
    x < MINEXP2(T) && return 0.0
 
    # reduce
    k = round(x) 
    n = _trunc(k)
    r = x - k

    # compute approximation
    x = _exp2(r)
    return _ldexp(x,n)
end