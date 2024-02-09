using LinearAlgebra

function get_tangent(A)
    # it is assumed that Jacobian A has full rank
    dim = size(A, 2) - size(A, 1)  # dimension of null space

    if dim > 0
        Q, _ = qr(A')
        z = Q[:, end - dim + 1:end]  # tangent space

        if det([A; z']) < 0 && dim == 1
            t = -z
        else
            t = z
        end
    elseif dim == 0
        t = []
    else
        error("Jacobian is of wrong dimension")
    end

    return t
end
