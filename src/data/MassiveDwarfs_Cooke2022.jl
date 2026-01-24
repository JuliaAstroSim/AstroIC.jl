"""
$(TYPEDSIGNATURES)

Load CO RC data for massive dwarfs from Cooke2022

Arguments:
- galaxy::String: Name of the dwarf galaxy

Returns:
- DataFrame with CO RC data including calculated velocity errors
"""
function load_massive_dwarf_CO_RC(galaxy::String)
    folder_data = @__DIR__
    df_CO_RC = DataFrame(CSV.File(joinpath(folder_data, "Cooke2022_RC", "$(galaxy).csv")))
    df_CO_RC.vel_e = df_CO_RC.vel_u .- df_CO_RC.vel
    df_CO_RC.vel_d = df_CO_RC.vel .- df_CO_RC.vel_e
    return df_CO_RC
end

"""
$(TYPEDSIGNATURES)

Load DM RC data for massive dwarfs from Cooke2022

Arguments:
- galaxy::String: Name of the dwarf galaxy

Returns:
- DataFrame with DM RC data including calculated velocity errors
"""
function load_massive_dwarf_DM_RC(galaxy::String)
    folder_data = @__DIR__
    df_DM_RC = DataFrame(CSV.File(joinpath(folder_data, "Cooke2022_RC_DM", "$(galaxy).csv")))
    df_DM_RC.vel_e = df_DM_RC.vel_u .- df_DM_RC.vel
    df_DM_RC.vel_d = df_DM_RC.vel .- df_DM_RC.vel_e
    return df_DM_RC
end
