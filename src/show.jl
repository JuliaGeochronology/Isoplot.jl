
# Custom pretty-printing for York fit results
function Base.show(io::IO, ::MIME"text/plain", x::YorkFit{T}) where T
    print(io, "YorkFit{$T}:\nLeast-squares linear fit of the form y = a + bx where")
    print(io, "\n  intercept a : $(x.intercept) (1σ)")
    print(io, "\n  slope b     : $(x.slope) (1σ)")
    print(io, "\n  MSWD        : $(x.mswd)\n")
end

Base.show(io::IO, ::MIME"text/plain", x::CI) = print(io, x)
function Base.print(io::IO, x::CI)
    l = round(x.mean - x.lower, sigdigits=3)
    u = round(x.upper - x.mean, sigdigits=3)
    d = floor(Int, log10(x.mean)) - floor(Int, log10(max(l,u)))
    m = round(x.mean, sigdigits=3+d)
    print(io, "$m +$u/-$l")
end
