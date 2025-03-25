@testset "broadcast                              " begin
    # example dimensions
    n = 10
    u = FTField(n, up_dealias_size(n))
    u .= rand.() .+im.*rand.()

    v = FTField(n, up_dealias_size(n))

    # return v and not the underlying array
    @test typeof(v .*= 1) <: FTField

    # use FTFields and scalar
    v .= u .+ conj.(u)
    for i in LinearIndices(v)
        @test imag(v[i]) == 0.0
    end

    v .= .- u .* u .* 2im
    for i in LinearIndices(v)
        @test v[i] == -u[i]*u[i]*2im
    end

    v .= 0.0 + 1im
    for i in LinearIndices(v)
        @test v[i] == 0.0 + 1im
    end

    v .-= 0.0 + 1im
    for i in LinearIndices(v)
        @test v[i] == 0.0 + 0.0*im
    end

    v .= 0.0
    v .+= u ./ u .+ u.^2
    for i in LinearIndices(v)
        @test v[i] == u[i]/u[i] + u[i]^2
    end
end

@testset "broadcast allocation                   " begin
    n = 10
    u = FTField(n, up_dealias_size(n))
    v = FTField(n, up_dealias_size(n))
    c = 1.0

    # use FTFields and scalar
    foo(v, u, c) = (v .= u .+ conj.(u) .* c)

    foo(v, u, c)
    @test (@allocated foo(v, u, c)) == 0
end