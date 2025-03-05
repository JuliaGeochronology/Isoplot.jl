
# Custom pretty-printing for York fit results
function Base.show(io::IO, x::YorkFit{T}) where T
    if get(io, :compact, false)
        print(io, "YorkFit{$T}(…)")
    else
        Base.show_default(io, x)
    end
end
function Base.show(io::IO, ::MIME"text/plain", x::YorkFit{T}) where T
    print(io, """YorkFit{$T}:
     Least-squares linear fit of the form y = a + bx with
      intercept: $(x.intercept) (1σ)
      slope    : $(x.slope) (1σ)
      MSWD     : $(x.mswd)
    """
    )
end

function Base.show(io::IO, x::CI{T}) where T
    if get(io, :compact, false)
        l = round(x.mean - x.lower, sigdigits=2)
        u = round(x.upper - x.mean, sigdigits=2)
        nodata = any(isnan, (x.mean, x.upper, x.lower))
        d = nodata ? 0 : floor(Int, log10(abs(x.mean))) - floor(Int, log10(max(abs(l),abs(u))))
        m = round(x.mean, sigdigits=3+d)
        print(io, "$m +$u/-$l")
    else
        Base.show_default(io, x)
    end
end
function Base.show(io::IO, ::MIME"text/plain", x::CI{T}) where T
    print(io, """CI{$T} $x (95% CI):
      mean  : $(x.mean)
      sigma : $(x.sigma)
      median: $(x.median)
      lower : $(x.lower)
      upper : $(x.upper)
    """
    )
end
function Base.print(io::IO, x::CI)
    l = round(x.mean - x.lower, sigdigits=3)
    u = round(x.upper - x.mean, sigdigits=3)
    nodata = any(isnan, (x.mean, x.upper, x.lower))
    d = nodata ? 0 : floor(Int, log10(abs(x.mean)) - log10(min(abs(l),abs(u))))
    m = round(x.mean, sigdigits=3+d)
    print(io, "$m +$u/-$l")
end
