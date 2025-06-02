# Overarching abstract type for raw data of all sorts derived from mass spectrometry
abstract type RawData{T<:AbstractFloat} <: Data{T} end 

# Concrete type for U-Pb SIMS data
struct UPbSIMSData{T<:AbstractFloat} <: RawData{T}
    Pb204::Vector{T}
    Pb206::Vector{T}
    Pb207::Vector{T}
    Pb208::Vector{T}
    Th232::Vector{T}
    U238::Vector{T}
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
end

"""
```julia
calibration(data::Collection{<:UPbSIMSData}, standardages::Collection; baseline::Number=0)
```
Create a `UPbSIMSCalibration` object `calib` given a dataset of 
standard SIMS analyses `data` with known ages `standardages`, to construct
a Pb/U RSF vs U/UO2 calibration line.

The resulting calibraiton line is stored in `calib.line`, while the raw data
is stored in `calib.data`.

### Examples
```
julia> standards = importsimsdata("./examples/data")
10-element Vector{Isoplot.UPbSIMSData{Float64}}:
 Isoplot.UPbSIMSData{Float64}(...)
 Isoplot.UPbSIMSData{Float64}(...)
 ⋮
 Isoplot.UPbSIMSData{Float64}(...)

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
function calibration(data::Collection{UPbSIMSData{T}}, standardages::Collection; baseline::Number=0) where {T<:AbstractFloat}
    standardratios = ratio.(standardages, λ238U)
    calib = similar(data, Analysis2D{T})
    for i in eachindex(data, standardratios)
        dᵢ = data[i]
        rUO2_U = (dᵢ.U238O2 .- baseline) ./ (dᵢ.U238 .- baseline)
        PbUrsf = (dᵢ.Pb206  .- baseline) ./ ((dᵢ.U238 .- baseline) .* value(standardratios[i]))
        calib[i] = Analysis(PbUrsf, rUO2_U)
    end
    return UPbSIMSCalibration(calib, yorkfit(calib))
end


"""
```julia
importsimsdata(dir::String; T=Float64)
```
Import all `.asc` files in the directory `dir` into a vector of `UPbSIMSData` objects.

### Examples
```
julia> standards = importsimsdata("./examples/data")
10-element Vector{Isoplot.UPbSIMSData{Float64}}:
 Isoplot.UPbSIMSData{Float64}(...)
 Isoplot.UPbSIMSData{Float64}(...)
 ⋮
 Isoplot.UPbSIMSData{Float64}(...)
```
"""
function importsimsdata(dir::String, kwargs...)
    @assert isdir(dir) "Expecting a directory"
    files = filter(x->contains(x, ".asc"), readdir(dir))
    @assert !isempty(files) "No .asc files found"
    return importsimsdata(joinpath.(dir, files), kwargs...)
end

function importsimsdata(files; T=Float64)
    @assert !isempty(files) "No files to import"
    data = similar(files, UPbSIMSData{T})
    for i in eachindex(files)
        data[i] = importsimsfile(files[i])
    end
    return data
end

function importsimsfile(filepath)
    @assert contains(filepath, ".asc") "Expecting .asc file"

    file = readdlm(filepath, '\t', skipblanks=true)
    datastart = findfirst(x->x=="RAW DATA:=======================================================================", file[:,1]) + 3
    dataend = findfirst(x->x=="PRIMARY INTENSITY DATA : ///////////////////////////////////////////////////////", file[:,1]) - 1
    hasdata = .!isempty.(file[datastart, :])

    # Check data labels
    datalabels = file[datastart-2, :]
    datalabels[3:end] = file[datastart-1, 2:end-1]
    datalabels = replace.(string.(datalabels), r"[ \t]+$" => "")
    @assert datalabels[5] == "204Pb"
    @assert datalabels[6] == "206Pb"
    @assert datalabels[7] == "207Pb"
    @assert datalabels[8] == "208Pb"
    @assert datalabels[9] == "232Th"
    @assert datalabels[10] == "238U"
    @assert datalabels[11] == "238U 16O2"

    data = Float64.(file[datastart:dataend, hasdata])
    return UPbSIMSData(data[:,5], data[:,6], data[:,7], data[:,8], data[:,9], data[:,10], data[:,11])
end


function age68(d::UPbSIMSData{T}, calib::UPbSIMSCalibration{T}; blank64::Number=0, baseline::Number=0) where {T}
    # Determine appropriate pb/u rsf correction
    rUO2_U = @. (d.U238O2 - baseline) / (d.U238 - baseline)
    PbUrsf = value(invline(calib.line, nanmean(rUO2_U)))
    
     r68 = @. ((d.Pb206 - baseline) - (d.Pb204 - baseline) * blank64) / ((d.U238 - baseline) * PbUrsf)
    map!(x-> x<0 ? T(NaN) : x, r68, r68)
    
    # Return 206Pb/238U age
    return @. log(1 + r68)/value(λ238U)
end
function age75(d::UPbSIMSData{T}, calib::UPbSIMSCalibration{T}; blank74::Number=0, U58::Number=1/137.818, baseline::Number=0) where {T}
    # Determine appropriate pb/u rsf correction
    rUO2_U = @. (d.U238O2 - baseline) / (d.U238 - baseline)
    PbUrsf = value(invline(calib.line, nanmean(rUO2_U)))
    
    # Calculate blank-, baseline-, and rsf-corrected 207/235 ratios
    r75 = @. ((d.Pb207 - baseline) - (d.Pb204 - baseline) * blank74) / ((d.U238 - baseline) * U58 * PbUrsf)
    map!(x-> x<0 ? T(NaN) : x, r75, r75)

    # Return 207Pb/235U age
    return @. log(1 +r75)/value(λ235U)
end


"""
```julia
calibrate(d::UPbSIMSData, calib::UPbSIMSCalibration, [cyclefilter]; 
    blank64::Number = first(stacey_kramers(0)), 
    blank74::Number = last(stacey_kramers(0)), 
    U58::Number = 1/137.818, 
    baseline::Number = 0,
)
```
Calibrate one or more U-Pb SIMS analyses `d` with the calibration `calib`, 
optionally filtering by `cyclefilter`, resulting in a `UPbAnalysis` for
each input `UPbSIMSData` object.

A Pb-204 based blank subtraction (common-Pb subtraction) is performed, 
using by default a present-day Stacey-Kramers common Pb composition.

U-235 is estimated based on measured U-238, assuming a 235/238 ratio 
of `U57` -- by default 1/137.818 (i.e., Heiss et al. 2012, doi: 10.1126/science.1215507)

A baseline-subtraction (equal for all masses) may optionally be performed 
given `baseline` excess counts per minute.

### Examples
```
julia> standards = importsimsdata("./examples/data")
10-element Vector{Isoplot.UPbSIMSData{Float64}}:
 Isoplot.UPbSIMSData{Float64}(...)
 Isoplot.UPbSIMSData{Float64}(...)
 ⋮
 Isoplot.UPbSIMSData{Float64}(...)

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
function calibrate(d::UPbSIMSData{T}, calib::UPbSIMSCalibration{T}, cf=:; blank64::Number=first(stacey_kramers(0)), blank74::Number=last(stacey_kramers(0)), U58::Number=1/137.818, baseline::Number=0) where {T}
    # Determine appropriate pb/u rsf correction
    rUO2_U = @. (d.U238O2 - baseline) / (d.U238 - baseline)
    PbUrsf = value(invline(calib.line, nanmean(rUO2_U)))
    
    # Calculate blank-, baseline-, and rsf-corrected 206/238 and 207/235 ratios
    r68 = @. ((d.Pb206 - baseline) - (d.Pb204 - baseline) * blank64) / ((d.U238 - baseline) * PbUrsf)
    r75 = @. ((d.Pb207 - baseline) - (d.Pb204 - baseline) * blank74) / ((d.U238 - baseline) * U58 * PbUrsf)

    # Return UPbAnalysis object
    return cf isa Colon ? UPbAnalysis(r75, r68) : UPbAnalysis(r75[cf], r68[cf])
end


"""
```julia
calibrate_blockwise(d::UPbSIMSData, calib::UPbSIMSCalibration, [cyclefilter]; 
    blank64::Number = first(stacey_kramers(0)), 
    blank74::Number = last(stacey_kramers(0)), 
    U58::Number = 1/137.818, 
    baseline::Number = 0,
    blocksize::Integer = 3,
)
```
Calibrate one or more U-Pb SIMS analyses `d` with the calibration `calib`, 
on a block-by-block basis into blocks of length `blocksize`, optionally 
filtering by `cyclefilter`, resulting in vector of `UPbAnalysis` objects
for each `UPbSIMSData` object.

A Pb-204 based blank subtraction (common-Pb subtraction) is performed, 
using by default a present-day Stacey-Kramers common Pb composition.

U-235 is estimated based on measured U-238, assuming a 235/238 ratio 
of `U57` -- by default 1/137.818 (i.e., Heiss et al. 2012, doi: 10.1126/science.1215507)

A baseline-subtraction (equal for all masses) may optionally be performed 
given `baseline` excess counts per minute.

### Examples
```
julia> standards = importsimsdata("./examples/data")
10-element Vector{Isoplot.UPbSIMSData{Float64}}:
 Isoplot.UPbSIMSData{Float64}(...)
 Isoplot.UPbSIMSData{Float64}(...)
 ⋮
 Isoplot.UPbSIMSData{Float64}(...)

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
function calibrate_blockwise(d::UPbSIMSData{T}, calib::UPbSIMSCalibration{T}; blank64::Number=first(stacey_kramers(0)), blank74::Number=last(stacey_kramers(0)), U58::Number=1/137.818, baseline::Number=0, blocksize::Integer=3) where {T}
    nblocks = length(d.U238)÷blocksize
    analyses = Vector{UPbAnalysis{T}}(undef, nblocks)
    for i in 1:nblocks
        bi = ((i-1)*blocksize+1):(i*blocksize)
        rUO2_U = @. (d.U238O2[bi] - baseline) / (d.U238[bi] - baseline)
        PbUrsf = value(invline(calib.line, nanmean(rUO2_U)))
        
        r68 = @. ((d.Pb206[bi] - baseline) - (d.Pb204[bi] - baseline) * blank64) / ((d.U238[bi] - baseline) * PbUrsf)
        r75 = @. ((d.Pb207[bi] - baseline) - (d.Pb204[bi] - baseline) * blank74) / ((d.U238[bi] - baseline) * U58 * PbUrsf)

        analyses[i] = UPbAnalysis(r75, r68)
    end
    return analyses
end
