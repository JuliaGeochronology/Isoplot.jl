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

struct UPbSIMSCalibration{T}
    data::Vector{BivariateAnalysis{T}}
    yf::YorkFit{T}
end

function calibrate(data::Collection{UPbSIMSData{T}}, standardages::Collection) where {T<:AbstractFloat}
    standardratios = ratio.(standardages, λ238U)
    calib = similar(data, BivariateAnalysis{T})
    for i in eachindex(data, standardratios)
        dᵢ = data[i]
        rUO2_U = dᵢ.U238O2 ./ dᵢ.U238
        PbUrsf = dᵢ.Pb206 ./ (dᵢ.U238 .* val(standardratios[i]))
        calib[i] = BivariateAnalysis(rUO2_U, PbUrsf)
    end
    return UPbSIMSCalibration(calib, yorkfit(calib))
end
export calibrate

function importsimsdata(dir)
    @assert isdir(dir) "Expecting a directory"
    files = filter(x->contains(x, ".asc"), readdir(dir))
    @assert !isempty(files) "No .asc files found"

    data = Vector{UPbSIMSData{Float64}}(undef, length(files))
    for i in eachindex(files)
        data[i] = importsimsfile(joinpath(dir, files[i]))
    end
    return data
end
export importsimsdata

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

# function reduce(d::UPbSIMSData, c::UPbSIMSCalibration)
# end