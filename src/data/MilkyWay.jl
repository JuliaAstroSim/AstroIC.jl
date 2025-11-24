"""
$(TYPEDSIGNATURES)
"""
function load_MW_RC_Eilers2019(filename = joinpath(@__DIR__, "MilkyWay/MW_RC_Eilers2019.dat"))
    dfEilers2019 = DataFrame(CSV.File(filename; header = false, delim=" ", ignorerepeated = true))
    rename!(dfEilers2019, [:r, :v, :σ_low, :σ_high])
    return dfEilers2019
end

"""
$(TYPEDSIGNATURES)
"""
function load_MW_RC_Mroz2019(filename = joinpath(@__DIR__, "MilkyWay/MW_RC_Mroz2019.dat"))
    dfMroz2019 = DataFrame(CSV.File(filename; header = false, delim=" ", ignorerepeated = true))
    rename!(dfMroz2019, [:r, :v, :σ_low, :σ_high])
    return dfMroz2019
end

"""
$(TYPEDSIGNATURES)
"""
function load_MW_RC_stddev_W21(filename = joinpath(@__DIR__, "MilkyWay/MW_RC_stddev_W21.txt"))
    dfstddev_W21 = DataFrame(CSV.File(filename; header = false, delim=" ", ignorerepeated = true))
    rename!(dfstddev_W21, [:id, :r, :v, :σ_v, :σ_r])
    return dfstddev_W21
end

"""
$(TYPEDSIGNATURES)
"""
function load_MW_RC_DS_W21(filename = joinpath(@__DIR__, "MilkyWay/MW_dsvel-W21-RC.txt"))
    dfDS_W21 = DataFrame(CSV.File(filename; header = false, skipto=3, delim="\t", ignorerepeated = true))
    rename!(dfDS_W21, [:r, :v_b, :v_CDM, :v_QUMOND, :v_MOG])
    return dfDS_W21
end