
denormalmin(::Type{Float64}) = 5.0e-324
denormalmin(::Type{Float32}) = 1.0f-45

# the following obviously rely on ^ being defined...
denormals(::Type{Float64}) = 2.0.^-(1023.0:1074.0) 
denormals(::Type{Float32}) = 2f0.^-(127f0:149f0)

tol = 5
@testset "Accuracy (in ulp) for $T" for T in (Float64, Float32)
    println("Accuracy tests for $T")
    
    xx = map(T, vcat(-10:0.0002:10, -1000:0.001:1000, -120:0.023:1000, -1000:0.02:2000))
    test_acc(T, Dict(Amal.exp => Base.exp), xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -120:0.023:1000, -1000:0.02:2000))
    test_acc(T, Dict(Amal.exp2 => Base.exp2), xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -35:0.023:1000, -300:0.01:300))
    test_acc(T, Dict(Amal.exp10 => Base.exp10), xx, tol)


    xx = map(T, vcat(0.0001:0.0001:10, 0.001:0.1:10000, 1.1.^(-1000:1000), 2.1.^(-1000:1000)))
    test_acc(T, Dict(Amal.log => Base.log), xx, tol)
    test_acc(T, Dict(Amal.log => Base.log), denormals(T), tol)


    # xx = map(T, vcat(0:0.2:10000, 1.1.^(-1000:1000), 2.1.^(-1000:1000)))
    # fun_table = Dict(Amal.cbrt => Base.cbrt)
    # test_acc(T, fun_table, xx, tol)

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
    # for i = 1:10
    #     s = reinterpret(T, reinterpret(inttype(T), T(pi)/4 * i) - asint(T,5000))
    #     e = reinterpret(T, reinterpret(inttype(T), T(pi)/4 * i) + asint(T,5000))
    #     d = s
    #     while d <= e 
    #         append!(xx, d)
    #         d = reinterpret(T, reinterpret(inttype(T), d) + asint(T,1))
    #     end
    # end
    # xx = append!(xx, -10:0.0002:10)
    # fun_table = Dict(Amal.tan => Base.tan)
    # test_acc(T, fun_table, xx, tol)

end
