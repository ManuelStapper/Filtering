function fitdims(H, R, F, Q, A = missing, z = missing)
    H = H .+ 0.0
    R = R .+ 0.0
    F = F .+ 0.0
    Q = Q .+ 0.0

    if typeof(R) == Float64
        R = reshape([R], (1, 1))
    end

    if !issymmetric(R)
        error("R must be a symmetric matrix.")
    else
        if det(R) <= 0
            error("R must be positive definite.")
        end
    end

    ny = size(R)[1]

    if typeof(Q) == Float64
        Q = reshape([Q], (1, 1))
    end

    if !issymmetric(Q)
        error("Q must be a symmetric matrix.")
    else
        if det(Q) <= 0
            error("Q must be positive definite.")
        end
    end

    nx = size(Q)[1]

    if typeof(H) == Vector{Float64}
        if length(H) != nx
            error("H must have suitable dimension.")
        end
        H = reshape(H, (1, nx))
    end


    if typeof(H) == Float64
        if nx != 1
            error("H must have suitable dimension.")
        end
        H = reshape([H], (1, 1))
    end

    if size(H) != (ny, nx)
        error("H must have suitable dimension.")
    end

    if typeof(F) == Float64
        if nx != 1
            error("H must have suitable dimension.")
        end
        F = reshape([F], (1, 1))
    end

    if size(F) != (nx, nx)
        error("F must have suitable dimension.")
    end

    if ismissing(z)
        if !ismissing(A)
            println("Specify A and z or none of them. Regressors are ignored.")
            A = missing
        end
    else
        if ismissing(A)
            println("Specify A and z or none of them. Regressors are ignored.")
            z = missing
        else
            A = A .+ 0.0
            z = z .+ 0.0

            if typeof(z) == Vector{Float64}
                z = reshape(z, (1, length(z)))
            end

            if typeof(z) != Matrix{Float64}
                error("z must be a suitable array.")
            end

            nz = size(z)[1]

            if typeof(A) == Vector{Float64}
                if (ny > 1) & (nz > 1)
                    error("A and z must have matching dimension.")
                end

                if length(A) == ny
                    A = reshape(A, (ny, 1))
                else
                    if length(A) == nz
                        A = reshape(A, (1, size(z)[1]))
                    else
                        error("A and z must have matching dimension.")
                    end
                end
            end

            if typeof(A) != Matrix{Float64}
                error("A must be a suitable array.")
            end

            if size(A) != (ny, nz)
                error("A must have suitable dimensions.")
            end
        end
    end

    return  H, R, F, Q, A, z
end
