# Simple linear regression

function linreg(x::AbstractVector{T}, y::AbstractVector{<:Number}) where {T<:Number}
    A = similar(x, length(x), 2)
    A[:,1] .= one(T)
    A[:,2] .= x
    return A\y
end


## -- The York (1968) two-dimensional linear regression with x and y uncertainties
    # as commonly used in isochrons

# Custom type to hold York fit resutls
struct YorkFit{T<:Number}
    intercept::Measurement{T}
    slope::Measurement{T}
    mswd::T
end

"""
```julia
yorkfit(x, σx, y, σy)
```
Uses the York (1968) two-dimensional least-squares fit to calculate `a`, `b`,
and uncertanties `σa`, `σb` for the equation `y = a + bx`, given `x`, `y` and
uncertaintes `σx`, ``σy`.

For further reference, see:
York, Derek (1968) "Least squares fitting of a straight line with correlated errors"
Earth and Planetary Science Letters 5, 320-324. doi: 10.1016/S0012-821X(68)80059-7

### Examples
```julia
julia> x = (1:100) .+ randn.();

julia> y = 2*(1:100) .+ randn.();

julia> yorkfit(x, ones(100), y, ones(100))
YorkFit{Float64}:
Least-squares linear fit of the form y = a + bx where
  intercept a : -0.29 ± 0.2 (1σ)
  slope b     : 2.0072 ± 0.0035 (1σ)
  MSWD        : 0.8136665223891004
```
"""
function yorkfit(x, σx, y, σy; niterations=10)

    ## 1. Ordinary linear regression (to get a first estimate of slope and intercept)

    # Check for missing data
    t = (x.==x) .& (y.==y) .& (σx.==σx) .& (σy.==σy)
    x = x[t]
    y = y[t]
    σx = σx[t]
    σy = σy[t]

    # Calculate the ordinary least-squares fit
    # For the equation y=a+bx, m(1)=a, m(2)=b
    a, b = linreg(x, y)

    ## 2. Now, let's define parameters needed by the York fit

    # Weighting factors
    ωx = 1.0 ./ σx.^2
    ωy = 1.0 ./ σy.^2

    # terms that don't depend on a or b
    α = sqrt.(ωx .* ωy)

    x̄ = sum(x)/length(x)
    ȳ = sum(y)/length(y)
    r = sum((x .- x̄).*(y .- ȳ)) ./ (sqrt(sum((x .- x̄).^2)) * sqrt(sum((y .- ȳ).^2)))

    ## 3. Perform the York fit (must iterate)
    W = ωx.*ωy ./ (b^2*ωy + ωx - 2*b*r.*α)

    X̄ = sum(W.*x) / sum(W)
    Ȳ = sum(W.*y) / sum(W)

    U = x .- X̄
    V = y .- Ȳ

    sV = W.^2 .* V .* (U./ωy + b.*V./ωx - r.*V./α)
    sU = W.^2 .* U .* (U./ωy + b.*V./ωx - b.*r.*U./α)
    b = sum(sV) ./ sum(sU)

    a = Ȳ - b .* X̄
    for i = 2:niterations
        W .= ωx.*ωy ./ (b^2*ωy + ωx - 2*b*r.*α)

        X̄ = sum(W.*x) / sum(W)
        Ȳ = sum(W.*y) / sum(W)

        U .= x .- X̄
        V .= y .- Ȳ

        sV .= W.^2 .* V .* (U./ωy + b.*V./ωx - r.*V./α)
        sU .= W.^2 .* U .* (U./ωy + b.*V./ωx - b.*r.*U./α)
        b = sum(sV) ./ sum(sU)

        a = Ȳ - b .* X̄
    end

    ## 4. Calculate uncertainties and MSWD
    β = W .* (U./ωy + b.*V./ωx - (b.*U+V).*r./α)

    u = X̄ .+ β
    v = Ȳ .+ b.*β

    xm = sum(W.*u)./sum(W)
    ym = sum(W.*v)./sum(W)

    σb = sqrt(1.0 ./ sum(W .* (u .- xm).^2))
    σa = sqrt(1.0 ./ sum(W) + xm.^2 .* σb.^2)

    # MSWD (reduced chi-squared) of the fit
    mswd = 1.0 ./ length(x) .* sum( (y .- a.-b.* x).^2 ./ (σy.^2 + b.^2 .* σx.^2) )

    ## Results
    return YorkFit(a ± σa, b ± σb, mswd)
end
