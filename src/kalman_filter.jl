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

        if typeof(x0) == Matrix{Float64}
            if size(x0)[1] != 1
                error("x0 needs to be a vector or 1-column matrix. Changed to default.")
                x0 = zeros(nx)
            else
                x0 = x0[:, 1]
            end
        end

        if typeof(x0) != Vector{Float64}
            error("Invalid input for x0. Changed to default.")
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
            if !isposdef(P0)
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
    xFilt = zeros(nx, T)
    PFilt = zeros(nx, nx, T)
    K = zeros(nx, ny, T)
    xPred = zeros(nx, T)
    PPred = zeros(nx, nx, T)

    xPred[:, 1] = F*x0
    PPred[:, :, 1] = F*P0*F' .+ Q

    K[:, :, 1] = (PPred[:, :, 1]*H')/(H*PPred[:, :, 1]*H' .+ R)
    if ismissing(z)
        xFilt[:, 1] = xPred[:, 1] + K[:, :, 1]*(y[:, 1] .- H*xPred[:, 1])
    else
        xFilt[:, 1] = xPred[:, 1] + K[:, :, 1]*(y[:, 1] .- A*z[:, 1] .- H*xPred[:, 1])
    end
    PFilt[:, :, 1] = PPred[:, :, 1] .- K[:, :, 1]*H*PPred[:, :, 1]

    for t = 2:T
        xPred[:, t] = F*xFilt[:, t-1]
        PPred[:, :, t] = F*PFilt[:, :, t-1]*F' .+ Q

        K[:, :, t] = (PPred[:, :, t]*H')/(H*PPred[:, :, t]*H' .+ R)
        if ismissing(z)
            xFilt[:, t] = xPred[:, t] + K[:, :, t]*(y[:, t] .- H*xPred[:, t])
        else
            xFilt[:, t] = xPred[:, t] + K[:, :, t]*(y[:, t] .- A*z[:, t] .- H*xPred[:, t])
        end
        PFilt[:, :, t] = PPred[:, :, t] .- K[:, :, t]*H*PPred[:, :, t]
    end

    return xFilt, PFilt
end
