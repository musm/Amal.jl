tol = 5

MRANGE(::Type{Float64}) = 10000000.0
MRANGE(::Type{Float32}) = 10000f0

inttype(::Type{Float64}) = Int64
inttype(::Type{Float32}) = Int32
inttype(::Type{Float16}) = Int16
asint{T}(::Type{T}, x) = inttype(T)(x)

@testset "Accuracy (in ulp) for $T" for T in (Float64, Float32 )
    println("Accuracy tests for $T")
    
    xx = map(T, vcat(-10:0.0002:10, -1000:0.001:1000, -120:0.023:1000, -1000:0.02:2000))
    fun_table = Dict(Amal.exp => Base.exp)
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -120:0.023:1000, -1000:0.02:2000))
    fun_table = Dict(Amal.exp2 => Base.exp2)
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -35:0.023:1000, -300:0.01:300))
    fun_table = Dict(Amal.exp10 => Base.exp10)
    test_acc(T, fun_table, xx, tol)


    # xx = map(T, vcat(-1000:0.021:1000, -1000:0.023:1000,  -10:0.0002:10,
    # -1:0.000002:1, 10.0.^-(0:0.02:300), -10.0.^-(0:0.02:300), 10.0.^(0:0.021:300), -10.0.^-(0:0.021:300)))
    # fun_table = Dict(Base.expm1 => Base.expm1)
    # test_acc(T, fun_table, xx, tol)

    # xx = map(T, vcat(-10:0.0002:10, -10000:0.2:10000, -10000:0.201:10000))
    # xx = map(T, vcat(-10:0.0002:10,-10:0.0001:10, -2:0.00001:2))
    # xx = map(T, vcat(-1:0.00001:1, -1:0.000001:1,-1:0.000021:1))
    # fun_table = Dict(Amal.atan => Base.atan)
    # test_acc(T, fun_table, xx, tol)


    # xx = T[]
    #     s = reinterpret(T, reinterpret(inttype(T), T(0)))
    #     e = reinterpret(T, reinterpret(inttype(T), T(pi)/4))
    #     d = s
    #     while d <= e 
    #         append!(xx, d)
    #         d = reinterpret(T, reinterpret(inttype(T), d) + asint(T,10))
    #     end
    # end
    # fun_table = Dict(Amal.tan => Base.tan)
    # test_acc(T, fun_table, xx, tol)


end
