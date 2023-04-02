
# Custom pretty-printing for York fit results
function Base.show(io::IO, ::MIME"text/plain", x::YorkFit{T}) where T
    print(io, """YorkFit{$T}:
     Least-squares linear fit of the form y = a + bx with
      intercept: $(x.intercept) (1σ)
      slope    : $(x.slope) (1σ)
      MSWD     : $(x.mswd)
    """
    )
end

function Base.show(io::IO, ::MIME"text/plain", x::CI{T}) where T
    print(io, """CI{$T} $x
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
    d = floor(Int, log10(abs(x.mean))) - floor(Int, log10(max(abs(l),abs(u))))
    m = round(x.mean, sigdigits=3+d)
    print(io, "$m +$u/-$l")
end
