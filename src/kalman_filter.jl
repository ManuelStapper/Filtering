function kalman_filter(y,
    H, R, F, Q,
    A = missing,
    z = missing,
    x0 = missing,
    P0 = missing)

    # Change input to arrays of floats
    H, R, F, Q, A, z = fitdims(H, R, F, Q, A, z)

    y = y .+ 0.0
    if typeof(y) == Vector{Float64}
        y = reshape(y, (1, length(y)))
    end

    if typeof(y) != Matrix{Float64}
        error("y not suitable.")
    end

    T = size(y)[2]
    ny = size(y)[1]
    nx = size(F)[1]
    if !ismissing(z)
        nz = size(z)[1]
        if size(z)[2] != T
            error("z must have T columns.")
        end
    else
        nz = missing
    end

    # Check initial distribution of states
    if ismissing(x0)
        x0 = zeros(nx)
    else
        x0 = x0 .+ 0.0
        if length(x0) == 1
            x0 = [x0]
        end

        if length(x0) != nx
            println("Dimension of x0 does not match nx, switched to x0 = 0.")
            x0 = zeros(nx)
        end
    end

    if !ismissing(P0)
        P0 = P0 .+ 0.0
        if length(P0) == 1
            if nx > 1
                println("Not suitable dimensions in P0, changed to default.")
                P0 = missing
            else
                if typeof(P0) == Float64
                    P0 = reshape([P0], (1, 1))
                end

                if typeof(P0) == Vector{Float64}
                    P0 = reshape(P0, (1, 1))
                end
            end
        end

        if !issymmetric(P0)
            println("P0 not symmetric, changed to default.")
            P0 = missing
        else
            if det(P0) <= 0
                println("P0 not positive definite, changed to default.")
                P0 = missing
            end
        end
    end

    if ismissing(P0)
        vecP0 = try
            (diagm(ones(nx^2)) .- kron(F, F))\vec(Q)
        catch
            println("Marginal variance for P0 failed, changed to unit variance.")
            vec(diagm(nx))
        end

        P0 = reshape(vecP0, (nx, nx))
    end

    # Start of Filtering recursions

end
