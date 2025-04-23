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

function calibration(data::Collection{UPbSIMSData{T}}, standardages::Collection) where {T<:AbstractFloat}
    standardratios = ratio.(standardages, λ238U)
    calib = similar(data, Analysis2D{T})
    for i in eachindex(data, standardratios)
        dᵢ = data[i]
        rUO2_U = dᵢ.U238O2 ./ dᵢ.U238
        PbUrsf = dᵢ.Pb206 ./ (dᵢ.U238 .* value(standardratios[i]))
        calib[i] = Analysis(PbUrsf, rUO2_U)
    end
    return UPbSIMSCalibration(calib, yorkfit(calib))
end

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


function age68(d::UPbSIMSData{T}, calib::UPbSIMSCalibration{T}; blank64::Number=0, baseline204::Number=0) where {T}
    # Determine appropriate pb/u rsf correction
    rUO2_U = d.U238O2 ./ d.U238
    PbUrsf = value(invline(calib.line, nanmean(rUO2_U)))
    
    # Calculate blank-, baseline-, and rsf-corrected 206/238 ratios
    r68 = @. (d.Pb206 - (d.Pb204 - baseline204) * blank64) / (d.U238 * PbUrsf)
    map!(x-> x<0 ? T(NaN) : x, r68, r68)
    
    # Return 206Pb/238U age
    return @. log(1 + r68)/value(λ238U)
end
function age75(d::UPbSIMSData{T}, calib::UPbSIMSCalibration{T}; blank74::Number=0, U58::Number=1/137.818, baseline204::Number=0) where {T}
    # Determine appropriate pb/u rsf correction
    rUO2_U = d.U238O2 ./ d.U238
    PbUrsf = value(invline(calib.line, nanmean(rUO2_U)))
    
    # Calculate blank-, baseline-, and rsf-corrected 207/235 ratios
    r75 = @. (d.Pb207 - (d.Pb204 - baseline204) * blank74) / (d.U238 * U58 * PbUrsf)
    map!(x-> x<0 ? T(NaN) : x, r75, r75)

    # Return 207Pb/235U age
    return @. log(1 +r75)/value(λ235U)
end



function calibrate(data::Collection{<:RawData}, calib::Calibration, cyclefilter=:; kwargs...)
    if (cyclefilter isa Collection) && (eachindex(cyclefilter) == eachindex(data))
        return [calibrate(data[i], calib, cyclefilter[i]; kwargs...) for i in eachindex(data)]
    else
        return [calibrate(data[i], calib, cyclefilter; kwargs...) for i in eachindex(data)]
    end
end
function calibrate(d::UPbSIMSData{T}, calib::UPbSIMSCalibration{T}, cf=:; blank64::Number=stacey_kramers(0)[1], blank74::Number=stacey_kramers(0)[2], U58::Number=1/137.818, baseline204::Number=0) where {T}
    # Determine appropriate pb/u rsf correction
    rUO2_U = d.U238O2 ./ d.U238
    PbUrsf = value(invline(calib.line, nanmean(rUO2_U)))
    
    # Calculate blank-, baseline-, and rsf-corrected 206/238 and 207/235 ratios
    r68 = @. (d.Pb206 - (d.Pb204 - baseline204) * blank64) / (d.U238 * PbUrsf)
    r75 = @. (d.Pb207 - (d.Pb204 - baseline204) * blank74) / (d.U238 * U58 * PbUrsf)

    # Return UPbAnalysis object
    return cf isa Colon ? UPbAnalysis(r75, r68) : UPbAnalysis(r75[cf], r68[cf])
end


function calibrate_blockwise(data::Collection{<:RawData}, calib::Calibration; kwargs...)
    return [calibrate_blockwise(data[i], calib; kwargs...) for i in eachindex(data)]
end
function calibrate_blockwise(d::UPbSIMSData{T}, calib::UPbSIMSCalibration{T}; blank64::Number=stacey_kramers(0)[1], blank74::Number=stacey_kramers(0)[2], U58::Number=1/137.818, baseline204::Number=0, blocksize=3) where {T}
    nblocks = length(d.U238)÷blocksize
    analyses = Vector{UPbAnalysis{T}}(undef, nblocks)
    for i in 1:nblocks
        bi = ((i-1)*blocksize+1):(i*blocksize)
        rUO2_U = d.U238O2[bi] ./ d.U238[bi]
        PbUrsf = value(invline(calib.line, nanmean(rUO2_U)))
        
        r68 = @. (d.Pb206[bi] - (d.Pb204[bi] - baseline204) * blank64) / (d.U238[bi] * PbUrsf)
        r75 = @. (d.Pb207[bi] - (d.Pb204[bi] - baseline204) * blank74) / (d.U238[bi] * U58 * PbUrsf)

        analyses[i] = UPbAnalysis(r75, r68)
    end
    return analyses
end
