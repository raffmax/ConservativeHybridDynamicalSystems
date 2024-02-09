using LinearAlgebra
using TensorOperations

function bifurcation_equation(u_BP, resFun, FiniteDifferenceStepSize)
    _, dH = resFun(u_BP)

    F, S, E = svd(dH, full = true)

    dim = 1

    E1 = E[:, end-1:end]
    E2 = E[:, 1:end-dim-1]
    F1 = F[:, end-dim+1:end]
    F2 = F[:, 1:end-dim]

    # central differences using TensorOperations
    ddH = zeros(size(dH, 1), size(dH, 2), dim+1)
    h = FiniteDifferenceStepSize

    @tensor begin
        for i = 1:dim+1
            ddH[:,:,i] = (resFun(u_BP + 0.5*h*E1[:, i])[2] - resFun(u_BP - 0.5*h*E1[:, i])[2]) / h
        end
    end

    Hessian = zeros(dim+1, dim+1)
    @tensor begin
        for i = 1:dim+1
            for j = 1:dim+1
                Hessian[i, j] = F1[:, 1]' * ddH[:, :, i] * E1[:, j]
            end
        end
    end

    # degenerate conic section
    d, V = eigen(Hessian)
    x1, y1 = 1, sqrt(-d[1] / d[2])
    x2, y2 = 1, -y1

    tang1 = [x1; y1]
    tang1 = V * tang1 / norm(tang1)
    tang2 = [x2; y2]
    tang2 = V * tang2 / norm(tang2)
    tangVec = hcat(E1*tang1, E1*tang2)

    return tangVec, Hessian
end