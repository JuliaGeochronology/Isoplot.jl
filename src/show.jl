
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
        minerr = min(abs(l),abs(u))
        nd = any(isnan, (x.mean, x.upper, x.lower)) || minerr == 0
        d = nd ? 0 : floor(Int, log10(abs(x.mean)) - log10(minerr))
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
    minerr = min(abs(l),abs(u))
    nd = any(isnan, (x.mean, x.upper, x.lower)) || minerr == 0
    d = nd ? 0 : floor(Int, log10(abs(x.mean)) - log10(minerr))
    m = round(x.mean, sigdigits=3+d)
    print(io, "$m +$u/-$l")
end

# Convert symbols to unicode
function prettify(s::Symbol)
    if s === :U238
        "²³⁸U"
    elseif s === :U238O
        "²³⁸U¹⁶O"
    elseif s === :U238O2
        "²³⁸U¹⁶O₂"
    elseif s === :U235
        "²³⁵U"
    elseif s === :Th232
        "²³²Th"
    elseif s === :Th232O
        "²³²Th¹⁶O"
    elseif s === :Th232O2
        "²³²Th¹⁶O₂"
    elseif s === :Pb208
        "²⁰⁸Pb"
    elseif s === :Pb207
        "²⁰⁷Pb"
    elseif s === :Pb206
        "²⁰⁶Pb"
    elseif s === :Pb204
        "²⁰⁴Pb"
    else
        # Any other element
        str = String(s)
        if contains(str, r"([A-Z][a-z]?)([0-9]*)")
            str = replace(str, r"([A-Z][a-z]?)([0-9]*)"=>s"\2\1")
        end
        superscriptnumbers(str)
    end
end

superscriptnumbers(s::AbstractString) = replace(s, "1"=>"¹", "2"=>"²", "3"=>"³", "4"=>"⁴", "5"=>"⁵", "6"=>"⁶", "7"=>"⁷", "8"=>"⁸", "9"=>"⁹", "0"=>"⁰", "+"=>"⁺", "-"=>"⁻",)
subscriptnumbers(s::AbstractString) = replace(s, "1"=>"₁", "2"=>"₂", "3"=>"₃", "4"=>"₄", "5"=>"₅", "6"=>"₆", "7"=>"₇", "8"=>"₈", "9"=>"₉", "0"=>"₀", "+"=>"₊", "-"=>"₋",)
