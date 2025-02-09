using DataInterpolations

function square(x::Real)
    return x^2
end

function smooth_data(x::AbstractVector, y::AbstractVector; spline_knots=5)
    fun = DataInterpolations.BSplineApprox(y, x, 3, spline_knots, :Uniform, :Uniform)
    return fun
end
