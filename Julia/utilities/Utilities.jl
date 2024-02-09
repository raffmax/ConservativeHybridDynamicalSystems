module Utilities

include("bifurcation_equation.jl")
export bifurcation_equation
#include("get_jacobian_fd.jl")
#export get_jacobian_fd
include("get_tangent.jl")
export get_tangent
include("newtons_method.jl")
export newtons_method
include("numerical_continuation.jl")
export numerical_continuation

end