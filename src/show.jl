
# Custom pretty-printing for York fit results
function Base.show(io::IO, ::MIME"text/plain", x::YorkFit{T}) where T
    print(io, "YorkFit{$T}:\nLeast-squares linear fit of the form y = a + bx where")
    print(io, "\n  intercept a : $(x.intercept) (1σ)")
    print(io, "\n  slope b     : $(x.slope) (1σ)")
    print(io, "\n  MSWD        : $(x.mswd)\n")
end
