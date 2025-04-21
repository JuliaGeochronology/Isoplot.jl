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

abstract type Calibration{T} end
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

function calibrate(data::Collection{<:RawData}, calib::Calibration; kwargs...)
    return [calibrate(data[i], calib; kwargs...) for i in eachindex(data)]
end
function calibrate(d::UPbSIMSData{T}, calib::UPbSIMSCalibration{T}; blank64::Number=stacey_kramers(0)[1], blank74::Number=stacey_kramers(0)[2], U58::Number=1/137.818, baseline204::Number=0) where {T}
    rUO2_U = d.U238O2 ./ d.U238
    PbUrsf = value(invline(calib.line, nanmean(rUO2_U)))
    
    r68 = @. (d.Pb206 - (d.Pb204 - baseline204) * blank64) / (d.U238 * PbUrsf)
    r75 = @. (d.Pb207 - (d.Pb204 - baseline204) * blank74) / (d.U238 * U58 * PbUrsf)

    return UPbAnalysis(r75, r68)
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