function kalman_sim(T::Int64, H, R, F, Q,
    A = missing,
    z = missing,
    burnin::Int64 = 500)

    H, R, F, Q, A, z = fitdims(H, R, F, Q, A, z)
    ny = size(H)[1]
    nx = size(F)[1]
    if !ismissing(z)
        nz = size(z)[1]
    else
        nz = missing
    end

    LR = cholesky(R).L
    ε = LR*rand(Normal(), (ny, T + burnin))

    LQ = cholesky(Q).L
    η = LQ*rand(Normal(), (nx, T + burnin))

    y = zeros(ny, T + burnin)
    x = zeros(nx, T + burnin)

    x[:, 1] = η[:, 1]
    y[:, 1] = H*x[:, 1] .+ ε[:, 1]

    for t = 2:burnin
        x[:, t] = F*x[:, t-1] .+ η[:, t]
        y[:, t] = H*x[:, t] .+ ε[:, t]
    end

    for t = burnin+1:burnin + T
        x[:, t] = F*x[:, t-1] .+ η[:, t]
        if ismissing(z)
            y[:, t] = H*x[:, t] .+ ε[:, t]
        else
            y[:, t] = A*z[:, t] .+ H*x[:, t] .+ ε[:, t]
        end
    end

    return x[:, burnin + 1:end], y[:, burnin + 1:end]
end
