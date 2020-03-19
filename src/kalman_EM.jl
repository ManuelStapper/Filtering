function kalman_EM(y, H, z = missing, x0 = missing, P0 = missing, critval::Float64 = 0.001)
    y = y .+ 0.0
    if typeof(y) == Vector{Float64}
        y = reshape(y, (1, length(y)))
    end

    if typeof(y) != Matrix{Float64}
        error("y needs to be a matrix.")
    end

    ny = size(y)[1]
    T = size(y)[2]

    H = H.+ 0.0

    if !ismissing(z)
        z = z + 0.0
        if typeof(z) == Vector{Float64}
            z = reshape(z, (1, T))
        end

        if typeof(z) != Matrix{Float64}
            error("z needs to be a matrix.")
        end

        nz = size(z)[1]

        if size(z)[2] != T
            error("Incorrect dimensions of z.")
        end
    end

    if typeof(H) == Float64
        if ny > 1
            error("Dimension of y and H not suitable.")
        end
        H = reshape([H], (1, 1))
    end

    if typeof(H) == Vector{Float64}
        if ny > 1
            if length(H) != ny
                error("Dimension of y and H not suitable.")
            end
            H = reshape(H, (ny, 1))
        else
            H = reshape(H, (1, length(H)))
        end
    end

    if typeof(H) != Matrix{Float64}
        error("Invalid input of H.")
    end

    nx = size(H)[2]

    if critval <= 0
        println("Invalid critical value, set to default.")
        critval = 0.001
    end

    if ismissing(x0)
        x0 = zeros(nx)
    else
        x0 = x0 .+ 0.0
        if length(x0) == 1
            x0 = [x0]
        end

        if length(x0) != nx
            println("Dimension of x0 does not match nx. x0 will be esimated.")
            x0 = missing
        end

        if typeof(x0) == Matrix{Float64}
            if size(x0)[1] != 1
                println("x0 needs to be a vector or 1-column matrix. x0 will be esimated")
                x0 = missing
            else
                x0 = x0[:, 1]
            end
        end

        if typeof(x0) != Vector{Float64}
            error("Invalid input for x0. x0 will be estimated.")
            x0 = missing
        end
    end

    if !ismissing(P0)
        P0 = P0 .+ 0.0
        if length(P0) == 1
            if nx > 1
                println("Not suitable dimensions in P0. P0 will be estimated.")
                P0 = missing
            end
        end

        if !issymmetric(P0)
            println("P0 not symmetric. P0 will be estimated")
            P0 = missing
        else
            if !isposdef(P0)
                println("P0 not positive definite. P0 will be estimated")
                P0 = missing
            end
        end
    end

    # If x0 or P0 missing, estimates are obtained from smoothing

    # Define initial values for matrices
    R = diagm(ones(ny))
    fvals = collect(1:nx)./(2*nx)
    F = diagm(fvals)
    Q = diagm(ones(nx))
    if !ismissing(z)
        A = fill(0.1, (ny, nz))
        Z = z[:, 1]*z[:, 1]'
        for t = 2:T
            Z .+= z[:, t]*z[:, t]
        end
        Z = inv(Z)
    else
        A = missing
    end



    # Smoothing results
    SR = kalman_smooth(y, H, R, F, Q, A, z, x0, P0)

    # Initialize EM-Steps
    checkval = Inf
    newQ = -1e9
    oldQ = -1e10
    while abs(newQ - oldQ) > critval
        xS = SR[3]
        PS = SR[4]
        PC = SR[5]
        x0S = SR[6]
        P0S = SR[7]

        AA = x0S*x0S' .+ P0S
        BB = xS[:, 1]*x0S' .+ PC[:, :, 1]
        CC = xS[:, 1]*xS[:, 1]' .+ PS[:, :, 1]

        for t = 2:T
            AA .+= xS[:, t-1]*xS[:, t-1]' .+ PS[:, :, t-1]
            BB .+= xS[:, t]*xS[:, t-1]' .+ PC[:, :, t]
            CC .+= xS[:, t]*xS[:, t]' .+ PS[:, :, t]
        end

        F = BB/AA
        Q = (CC .- F*BB')./T
        if !ismissing(A)
            A = sum((y .- H*xS), dims = 2)'*Z
        end

        u = y .- H*xS
        if !ismissing(z)
            u .+= A*z
        end

        temp = u[:, 1]*u[:, 1]' .+ H*PS[:, :, 1]*H'
        for t = 2:T
            temp .+= u[:, t]*u[:, t]' .+ H*PS[:, :, t]*H'
        end
        R = temp./T

        oldQ = copy(newQ)
        newQ = -T/2*log(det(Q)) - T/2*log(det(R))
        newQ -= 1/2*tr(inv(Q)*(CC .- BB*F' .- F*BB' + F*AA*F'))
        newQ -= 1/2*tr(inv(R)*temp)

        SR = kalman_smooth(y, H, R, F, Q, A, z, ifelse(ismissing(x0), x0S, x0), ifelse(ismissing(P0), P0S, P0))
    end

    return H, R, F, Q, A, ifelse(ismissing(x0), x0S, x0), ifelse(ismissing(P0), P0S, P0)
end
