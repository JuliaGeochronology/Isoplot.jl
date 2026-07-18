# Overarching abstract type for raw data of all sorts derived from mass spectrometry
abstract type RawData{T<:AbstractFloat} <: Data{T} end 

# Concrete type for U-Pb SIMS data
struct UThPbSIMSData{T<:AbstractFloat} <: RawData{T}
    Pb204::Vector{T}
    Pb206::Vector{T}
    Pb207::Vector{T}
    Pb208::Vector{T}
    Th232::Vector{T}
    Th232O::Vector{T}
    Th232O2::Vector{T}
    U238::Vector{T}
    U238O::Vector{T}
    U238O2::Vector{T}
end

# Overarching abstract type for calibrations of all sorts
abstract type Calibration{T} end

# Generic methods to allow broadcasting and comparison
Base.length(x::Calibration) = 1
Base.iterate(x::Calibration) = (x, nothing)
Base.iterate(x::Calibration, state) = nothing
Base.:(==)(x::Calibration, y::Calibration) = false
function Base.:(==)(x::T, y::T) where {T<:Calibration}
    for n in fieldnames(T)
        isequal(getfield(x, n), getfield(y, n)) || return false
    end
    return true
end

struct UPbSIMSCalibration{T} <: Calibration{T}
    data::Vector{Analysis2D{T}}
    line::YorkFit{T}
    numerator::Symbol
    denominator::Symbol
    U::Symbol
    Th::Symbol
end

"""
```julia
calibration(data::Collection{UThPbSIMSData{T}}, standardages::Collection; 
    \tblank = stacey_kramers,       # Common Pb composition in order (206/204, 207/204, 208/204). If the stacey_kramers function itself is given without argument, the age of each standard will be used.
    \tcorrection = :none,           # [:none, :Pb204, :Pb208]
    \titerations = 8,               # Number of iterations for Pb-208 correction
    \tThUrsf::Number = 1.0,         # For Pb-208 correction; default is no fractionation
)
```
Create a `UPbSIMSCalibration` object `calib` given a dataset of 
standard SIMS analyses `data` with known ages `standardages`, to construct
a Pb/U RSF vs U/UO2 calibration line.

The resulting calibraiton line is stored in `calib.line`, while the raw data
is stored in `calib.data`.

### Examples
```
julia> standards = importsimsdata("./examples/data")
10-element Vector{Isoplot.UThPbSIMSData{Float64}}:
 Isoplot.UThPbSIMSData{Float64}(...)
 Isoplot.UThPbSIMSData{Float64}(...)
 ⋮
 Isoplot.UThPbSIMSData{Float64}(...)

julia> standardages = fill(1099., length(standards));

julia> calib = calibration(standards, standardages)
Isoplot.UPbSIMSCalibration{Float64}(...)

julia> calib.line
YorkFit{Float64}:
 Least-squares linear fit of the form y = a + bx with
  intercept: 1.25 ± 0.19 (1σ)
  slope    : 0.82 ± 0.066 (1σ)
  MSWD     : 0.2704923359277713
```
"""
function calibration(data::Collection{UThPbSIMSData{T}}, standardages::Collection; 
        blank::Union{NTuple{3,Number}, Function} = stacey_kramers,
        correction::Symbol = :none,
        iterations::Integer = 8,
        numerator::Symbol = :U238O2,
        denominator::Symbol = :U238,
        U::Symbol = :U238,
        Th::Symbol = :Th232,
    ) where {T<:AbstractFloat}
    @assert correction ∈ (:none, :Pb204, :Pb208,) "Common Pb correction options are `:none`, `:Pb204`, or `:Pb208`"
    @assert U ∈ (:U238, :U238O, :U238O2,)
    @assert Th ∈ (:Th232, :Th232O, :Th232O2,)
    # Initialize blank ratios
    blank64, blank74, blank84 = if blank isa NTuple
        blank
    else
        0., 0., 0.
    end

    # Allocate array for analyses
    calib = similar(data, Analysis2D{T})

    if correction === :Pb204
        # Pb-204 corection
        for i in eachindex(data)
            dᵢ = data[i]
            standardratioPb206U238 = ratio(standardages[i], value(λ238U))
            if blank isa Function
                blank64, blank74, blank84 = blank(standardages[i])
            end
            oxideratio = getfield(dᵢ, numerator) ./ getfield(dᵢ, denominator)
            Pb206c = dᵢ.Pb204 .* blank64
            PbUrsf = (dᵢ.Pb206 .- Pb206c) ./ (getfield(dᵢ, U) .* standardratioPb206U238)
            calib[i] = Analysis(PbUrsf, oxideratio)
        end
        
    elseif correction === :Pb208
        # Pb-208 corection
        for i in eachindex(data)
            dᵢ = data[i]
            standardratioPb206U238 = ratio(standardages[i], value(λ238U))
            standardratioPb208Th232 = ratio(standardages[i], value(λ232Th))
            if blank isa Function
                blank64, blank74, blank84 = blank(standardages[i])
            end
            blank68 = blank64/blank84
            U238 = getfield(dᵢ, U)
            Th232 = getfield(dᵢ, Th)
            oxideratio = getfield(dᵢ, numerator) ./ getfield(dᵢ, denominator)
            Pb208r, Pb206c = zeros(size(oxideratio)), zeros(size(oxideratio))
            PbUrsf = dᵢ.Pb206 ./ (getfield(dᵢ, U) .* standardratioPb206U238) # No subtraction the first time
            for _ in 1:iterations
                @. Pb208r = Th232 * standardratioPb208Th232 * PbUrsf
                @. Pb206c = (dᵢ.Pb208 - Pb208r) * blank68
                @. PbUrsf = (dᵢ.Pb206 - Pb206c) / (U238 * standardratioPb206U238)
            end
            calib[i] = Analysis(PbUrsf, oxideratio)
        end

    else # :none
        # No common Pb correction
        for i in eachindex(data)
            dᵢ = data[i]
            standardratioPb206U238 = ratio(standardages[i], value(λ238U))
            oxideratio = getfield(dᵢ, numerator) ./ getfield(dᵢ, denominator)
            PbUrsf = dᵢ.Pb206 ./ (getfield(dᵢ, U) .* standardratioPb206U238)
            calib[i] = Analysis(PbUrsf, oxideratio)
        end
    end
    return UPbSIMSCalibration(calib, yorkfit(calib), numerator, denominator, U, Th)
end


"""
```julia
importsimsdata(dir::String; T=Float64, baseline=0)
```
Import all `.asc` files in the directory `dir` into a vector of `UThPbSIMSData` objects.

A baseline-subtraction (equal for all masses) may optionally be performed 
given `baseline` excess counts per minute.

### Examples
```
julia> standards = importsimsdata("./examples/data")
10-element Vector{Isoplot.UThPbSIMSData{Float64}}:
 Isoplot.UThPbSIMSData{Float64}(...)
 Isoplot.UThPbSIMSData{Float64}(...)
 ⋮
 Isoplot.UThPbSIMSData{Float64}(...)
```
"""
function importsimsdata(dir::String; kwargs...)
    @assert isdir(dir) "Expecting a directory"
    files = filter(x->contains(x, ".asc"), readdir(dir))
    @assert !isempty(files) "No .asc files found"
    return importsimsdata(joinpath.(dir, files); kwargs...)
end
function importsimsdata(files::AbstractVector; T=Float64, kwargs...)
    @assert !isempty(files) "No files to import"
    data = Vector{UThPbSIMSData{T}}(undef, length(files))
    for i in eachindex(data)
        data[i] = importsimsfile(files[i]; kwargs...)
    end
    return data
end

function importsimsfile(filepath::AbstractString; 
        system::Symbol=:UThPb,
        baseline::Number=0,
    )
    @assert system ∈ (:UThPb,) "System $system not recognized"
    @assert contains(filepath, ".asc") "Expecting .asc file"

    # Read data file
    file = readdlm(filepath, '\t', skipblanks=true)
    datastartmarker = findfirst(x->x=="RAW DATA:=======================================================================", file[:,1])
    @assert !isnothing(datastartmarker) "Could not find data start in $filepath"
    datastart = datastartmarker + 3
    dataendmarker = findfirst(x->x=="PRIMARY INTENSITY DATA : ///////////////////////////////////////////////////////", file[:,1])
    @assert !isnothing(dataendmarker) "Could not find data start in $filepath"
    dataend = dataendmarker - 1
    @assert datastart <= dataend "No data rows found in $filepath"
    hasdata = .!isempty.(file[datastart, :])
    @assert any(hasdata) "No data columns foun in $filepath"
    data = Float64.(file[datastart:dataend, hasdata])

    # Check data labels
    datalabels = file[datastart-2, :]
    datalabels[3:end] = file[datastart-1, 2:end-1]
    datalabels = replace.(string.(datalabels[hasdata]), r"[ \t]+$" => "")

    # Pb isotopes
    Pb204 = findcounts("204Pb", datalabels, data)
    Pb206 = findcounts("206Pb", datalabels, data)
    Pb207 = findcounts("207Pb", datalabels, data)
    Pb208 = findcounts("208Pb", datalabels, data)
    # Th and Th oxides
    Th232 = findcounts("232Th", datalabels, data)
    Th232O = findcounts("232Th 16O", datalabels, data)
    Th232O2 = findcounts("232Th 16O2", datalabels, data)
    # U and U oxides
    U238 = findcounts("238U", datalabels, data)
    U238O = findcounts("238U 16O", datalabels, data)
    U238O2 = findcounts("238U 16O2", datalabels, data)

    # Baseline subtraction, if any
    if baseline != 0
        # Smallest observed count rate
        minobserved = nanminimum(nanmean.((Pb204, Pb206, Pb207, Pb208, Th232, Th232O, Th232O2, U238, U238O, U238O2)))
        # Subtract, ensuring that mean count rate cannot ever become negative
        maxbaseline = min(baseline, minobserved)
        (Pb204, Pb206, Pb207, Pb208, Th232, Th232O, Th232O2, U238, U238O, U238O2) .|> v -> v.-=maxbaseline
    end

    return UThPbSIMSData(Pb204, Pb206, Pb207, Pb208, Th232, Th232O, Th232O2, U238, U238O, U238O2)
end
function findcounts(label, datalabels, data::AbstractMatrix{T}, missingval::T=T(NaN)) where {T<:Number}
    if label ∈ datalabels
        col = findfirst(isequal(label), datalabels)
        return data[:,col]
    else
        return fill(missingval, size(data, 1))
    end
end


function age68(d::UThPbSIMSData{T}, calib::UPbSIMSCalibration{T}; blank64::Number=0) where {T}
    # Determine appropriate pb/u rsf correction
    oxideratio = getfield(d, calib.numerator) ./ getfield(d, calib.denominator)
    PbUrsf = value(invline(calib.line, nanmean(oxideratio)))
    
    r68 = @. (d.Pb206 - d.Pb204 * blank64) / (d.U238 * PbUrsf)
    map!(x-> x<0 ? T(NaN) : x, r68, r68)
    
    # Return 206Pb/238U age
    return @. log(1 + r68)/value(λ238U)
end
function age75(d::UThPbSIMSData{T}, calib::UPbSIMSCalibration{T}; blank74::Number=0, U58::Number=1/137.818) where {T}
    # Determine appropriate pb/u rsf correction
    oxideratio = getfield(d, calib.numerator) ./ getfield(d, calib.denominator)
    PbUrsf = value(invline(calib.line, nanmean(oxideratio)))
    
    # Calculate blank-, and rsf-corrected 207/235 ratios
    r75 = @. (d.Pb207 - d.Pb204 * blank74) / (d.U238 * U58 * PbUrsf)
    map!(x-> x<0 ? T(NaN) : x, r75, r75)

    # Return 207Pb/235U age
    return @. log(1 +r75)/value(λ235U)
end


"""
```julia
calibrate(d::UThPbSIMSData, calib::UPbSIMSCalibration, [cyclefilter]; 
    \tblank = stacey_kramers(0),
    \tU58 = 1/137.818, 
    \tcorrection = :Pb204,
    \titerations = 8
)
```
Calibrate one or more U-Pb SIMS analyses `d` with the calibration `calib`, 
optionally filtering by `cyclefilter`, resulting in a `UPbAnalysis` for
each input `UThPbSIMSData` object.

A Pb-204 based blank subtraction (common-Pb subtraction) is performed, 
using by default a present-day Stacey-Kramers common Pb composition.

U-235 is estimated based on measured U-238, assuming a 235/238 ratio 
of `U58` -- by default 1/137.818 (i.e., Heiss et al. 2012, doi: 10.1126/science.1215507)

### Examples
```
julia> standards = importsimsdata("./examples/data")
10-element Vector{Isoplot.UThPbSIMSData{Float64}}:
 Isoplot.UThPbSIMSData{Float64}(...)
 Isoplot.UThPbSIMSData{Float64}(...)
 ⋮
 Isoplot.UThPbSIMSData{Float64}(...)

julia> standardages = fill(1099., length(standards));

julia> calib = calibration(standards, standardages)
Isoplot.UPbSIMSCalibration{Float64}(...)

julia> calibrate(standards[1], calib)
UPbAnalysis{Float64}(Analysis2D{Float64}([1.8930052834566204, 0.1835092671444941], [0.021732106262365054, 0.0011423135956512402], [0.0004722844425987264 1.4604632531997064e-5; 1.4604632531997064e-5 1.304880350809665e-6]))
```
"""
function calibrate(data::Collection{<:RawData}, calib::Calibration, cyclefilter=:; kwargs...)
    if (cyclefilter isa Collection) && (eachindex(cyclefilter) == eachindex(data))
        return [calibrate(data[i], calib, cyclefilter[i]; kwargs...) for i in eachindex(data)]
    else
        return [calibrate(data[i], calib, cyclefilter; kwargs...) for i in eachindex(data)]
    end
end
function calibrate(d::UThPbSIMSData{T}, calib::UPbSIMSCalibration{T}, cf=:; 
        blank::NTuple{3,Number} = stacey_kramers(0),
        U58::Number = 1/137.818, 
        correction::Symbol = :Pb204,
        iterations::Integer = 8,
    ) where {T}
    @assert correction ∈ (:none, :Pb204, :Pb208) "Common Pb correction options are `:none`, `:Pb204`, or `:Pb208`"

    # Determine appropriate pb/u rsf correction
    oxideratio = getfield(d, calib.numerator) ./ getfield(d, calib.denominator)
    PbUrsf = value(invline(calib.line, nanmean(oxideratio)))
    U238 = getfield(d, calib.U)
    Th232 = getfield(d, calib.Th)

    # Calculate blank- and rsf-corrected 206/238 and 207/235 ratios
    blank64, blank74, blank84 = blank
    if correction === :Pb204
        # Pb-204 correction
        r68 = @. (d.Pb206 - d.Pb204 * blank64) / (U238 * PbUrsf)
        r75 = @. (d.Pb207 - d.Pb204 * blank74) / (U238 * U58 * PbUrsf)
        r82 = @. (d.Pb208 - d.Pb204 * blank84) / (Th232 * PbUrsf)

    elseif correction === :Pb208
        # Non-iterative Pb-208 correction
        blank68, blank78 = blank64/blank84, blank74/blank84
        r68 = @. (d.Pb206 - d.Pb208 * blank68) / (U238 * PbUrsf)
        r75 = @. (d.Pb207 - d.Pb208 * blank78) / (U238 * U58 * PbUrsf)
        r82 = @. d.Pb208 * NaN
        Pb208r, agest = zeros(size(r68)), zeros(size(r68))
        for _ in 1:iterations
            @. agest = log(1 + max(r68, zero(T)))/value(λ238U)
            @. Pb208r = Th232 * ratio(agest, value(λ232Th)) * PbUrsf
            @. r68 = (d.Pb206 - (d.Pb208 - Pb208r) * blank68) / (U238 * PbUrsf)
            @. r75 = (d.Pb207 - (d.Pb208 - Pb208r) * blank78) / (U238 * U58 * PbUrsf)
        end

    else # none
        # No common Pb correction
        r68 = @. d.Pb206 / (U238 * PbUrsf)
        r75 = @. d.Pb207 / (U238 * U58 * PbUrsf)
        r82 = @. d.Pb208 / (Th232 * PbUrsf)

    end
    
    # Return UPbAnalysis object
    return cf isa Colon ? UThPbAnalysis(r75, r68, r82) : UThPbAnalysis(r75[cf], r68[cf], r82[cf])
end


"""
```julia
calibrate_blockwise(d::UThPbSIMSData, calib::UPbSIMSCalibration, [cyclefilter]; 
    \tblocksize::Integer = 3,
    \tblank::NTuple{3,Number} = stacey_kramers(0),
    \tU58::Number = 1/137.818, 
    \tcorrection::Symbol = :Pb204,
    \titerations::Integer = 8,
)
```
Calibrate one or more U-Pb SIMS analyses `d` with the calibration `calib`, 
on a block-by-block basis into blocks of length `blocksize`, optionally 
filtering by `cyclefilter`, resulting in vector of `UPbAnalysis` objects
for each `UThPbSIMSData` object.

A Pb-204 based blank subtraction (common-Pb subtraction) is performed, 
using by default a present-day Stacey-Kramers common Pb composition.

U-235 is estimated based on measured U-238, assuming a 235/238 ratio 
of `U57` -- by default 1/137.818 (i.e., Heiss et al. 2012, doi: 10.1126/science.1215507)

### Examples
```
julia> standards = importsimsdata("./examples/data")
10-element Vector{Isoplot.UThPbSIMSData{Float64}}:
 Isoplot.UThPbSIMSData{Float64}(...)
 Isoplot.UThPbSIMSData{Float64}(...)
 ⋮
 Isoplot.UThPbSIMSData{Float64}(...)

julia> standardages = fill(1099., length(standards));

julia> calib = calibration(standards, standardages)
Isoplot.UPbSIMSCalibration{Float64}(...)

julia> calibrate_blockwise(standards[1], calib)
6-element Vector{UPbAnalysis{Float64}}:
 UPbAnalysis{Float64}(Analysis2D{Float64}([1.9522725644621401, 0.19479933826738347], [0.10870961850778198, 0.0033669214494005937], [0.011817781156107494 0.0003491283425207373; 0.0003491283425207373 1.1336160046433794e-5]))
 UPbAnalysis{Float64}(Analysis2D{Float64}([1.970015175457096, 0.18617712628434524], [0.01431116289886847, 0.0023768591910531384], [0.00020480938351794938 -1.1409841282689287e-6; -1.1409841282689287e-6 5.649459614093779e-6]))
 UPbAnalysis{Float64}(Analysis2D{Float64}([1.868382104433831, 0.18258981150055945], [0.03784079639497292, 0.0012781398026132764], [0.0014319258718057954 -4.8258045726099764e-5; -4.8258045726099764e-5 1.633641355024305e-6]))
 UPbAnalysis{Float64}(Analysis2D{Float64}([1.8869190705183156, 0.18091237376063943], [0.038037522681535894, 0.0013370053746898942], [0.0014468531317483576 -1.2603028339239038e-6; -1.2603028339239038e-6 1.7875833719496643e-6]))
 UPbAnalysis{Float64}(Analysis2D{Float64}([1.856977011176083, 0.17882174298755873], [0.07424037726958133, 0.0014834460179062856], [0.005511633617129768 8.764295590103794e-5; 8.764295590103794e-5 2.2006120880420158e-6]))
 UPbAnalysis{Float64}(Analysis2D{Float64}([1.8386867975254575, 0.18063423529792652], [0.05211549498674472, 0.00248005466136501], [0.002716024817713414 9.982927603459965e-5; 9.982927603459965e-5 6.150671123358313e-6]))
```
"""
function calibrate_blockwise(data::Collection{<:RawData}, calib::Calibration; kwargs...)
    return [calibrate_blockwise(data[i], calib; kwargs...) for i in eachindex(data)]
end
function calibrate_blockwise(d::UThPbSIMSData{T}, calib::UPbSIMSCalibration{T}; 
        blocksize::Integer = 3,
        blank::NTuple{3,Number} = stacey_kramers(0),
        U58::Number = 1/137.818, 
        correction::Symbol = :Pb204,
        iterations::Integer = 8,
    ) where {T}
    @assert correction ∈ (:none, :Pb204, :Pb208) "Common Pb correction options are `:none`, `:Pb204`, or `:Pb208`"
    blank64, blank74, blank84 = blank

    # Allocate result array
    nblocks = length(d.U238)÷blocksize
    analyses = Vector{UThPbAnalysis{T}}(undef, nblocks)

    # Cycle through each block
    for i in 1:nblocks
        # Determine appropriate pb/u rsf correction for this block
        bi = ((i-1)*blocksize+1):(i*blocksize)
        oxideratio = getfield(d, calib.numerator)[bi] ./ getfield(d, calib.denominator)[bi]
        PbUrsf = value(invline(calib.line, nanmean(oxideratio)))
        U238 = getfield(d, calib.U)
        Th232 = getfield(d, calib.Th)

        # Calculate blank- and rsf-corrected 206/238 and 207/235 ratios
        if correction === :Pb204
            # Pb-204 correction
            r68 = @. (d.Pb206[bi] - d.Pb204[bi] * blank64) / (U238[bi] * PbUrsf)
            r75 = @. (d.Pb207[bi] - d.Pb204[bi] * blank74) / (U238[bi] * U58 * PbUrsf)
            r82 = @. (d.Pb208[bi] - d.Pb204[bi] * blank84) / (Th232[bi] * PbUrsf)

        elseif correction === :Pb208
            # Non-iterative Pb-208 correction
            blank68, blank78 = blank64/blank84, blank74/blank84
            r68 = @. (d.Pb206[bi] - d.Pb208[bi] * blank68) / (U238[bi] * PbUrsf)
            r75 = @. (d.Pb207[bi] - d.Pb208[bi] * blank78) / (U238[bi] * U58 * PbUrsf)
            r82 = fill(NaN, length(bi))
            Pb208r, agest = zeros(size(r68)), zeros(size(r68))
            for _ in 1:iterations
                @. agest = log(1 + max(r68, zero(T)))/value(λ238U)
                @. Pb208r = Th232[bi] * ratio(agest, value(λ232Th)) * PbUrsf
                @. r68 = (d.Pb206[bi] - (d.Pb208[bi] - Pb208r) * blank68) / (U238[bi] * PbUrsf)
                @. r75 = (d.Pb207[bi] - (d.Pb208[bi] - Pb208r) * blank78) / (U238[bi] * U58 * PbUrsf)
            end

        else # none
            # No common Pb correction
            r68 = @. d.Pb206[bi] / (U238[bi] * PbUrsf)
            r75 = @. d.Pb207[bi] / (U238[bi] * U58 * PbUrsf)
            r82 = @. d.Pb208[bi] / (Th232[bi] * PbUrsf)

        end
        
        analyses[i] = UThPbAnalysis(r75, r68, r82)
    end
    return analyses
end
