function newtons_method(fun, x, options; varargin...)
    if length(varargin) > 0
        f, jacobian, t = varargin
        fun_counter = 0
    else
        f, jacobian = fun(x)
        fun_counter = 1
        t = get_tangent(jacobian)
    end

    err = maximum(abs.(f))
    iter = 0

    while err > options[:FunctionTolerance] && iter < options[:MaxIterations]
        if err > 1
            break
        end

        if options[:aimOnTarget]
            idx = options[:idxConPar]
            delta_x = -[jacobian; hcat(zeros(1, idx - 1), 1, zeros(1, length(x) - idx))] \ [f; 0]
        else
            delta_x = -[jacobian; t'] \ [f; zeros(size(t, 2))]
        end

        if norm(delta_x, 1) < options[:StepTolerance]
            println("The Newton step is smaller than StepTolerance ", options[:StepTolerance], "!")
            break
        end

        x = delta_x + x

        f, jacobian = fun(x)
        fun_counter += 1

        t = get_tangent(jacobian)
        iter += 1
        err = maximum(abs.(f))
    end

    fval = f
    if err < options[:FunctionTolerance]
        exitflag = 1
    else
        exitflag = -1
    end

    output = Dict(
        :t => get_tangent(jacobian),
        :iterations => iter,
        :funcCount => fun_counter
    )

    return x, fval, exitflag, output, jacobian
end