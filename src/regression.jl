# Simple linear regression

function linreg(x::AbstractVector{T}, y::AbstractVector{<:Number}) where {T<:Number}
    A = similar(x, length(x), 2)
    A[:,1] .= one(T)
    A[:,2] .= x
    return A\y
end
