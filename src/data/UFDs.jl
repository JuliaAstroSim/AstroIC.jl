"""
$(TYPEDSIGNATURES)

Load `src/data/UFDs.csv` with units

Reference: Table 1 and Table A.1 of hayashi2023dark_Astrophys. J. - Dark matter halo properties of the galactic dwarf satellites implication for chemo-dynamical evolution
"""
function load_data_UFDs()
    df_UFDs = DataFrame(CSV.File(joinpath(@__DIR__, "UFDs.csv"), delim=" ", skipto=3, header=2, ignorerepeated = true));
    df_UFDs.L   = (10 .^ df_UFDs.log10L) * u"Lsun"
    df_UFDs.L_u = (10 .^ (df_UFDs.log10L + df_UFDs.log10L_u)) * u"Lsun"
    df_UFDs.L_d = (10 .^ (df_UFDs.log10L + df_UFDs.log10L_d)) * u"Lsun"
    df_UFDs.b   = (10 .^ df_UFDs.log10b) * u"pc"
    df_UFDs.b_u = (10 .^ (df_UFDs.log10b + df_UFDs.log10b_u)) * u"pc"
    df_UFDs.b_d = (10 .^ (df_UFDs.log10b + df_UFDs.log10b_d)) * u"pc"
    df_UFDs.rho0   = (10 .^ df_UFDs.log10rho0) * u"Msun/pc^3"
    df_UFDs.rho0_u = (10 .^ (df_UFDs.log10rho0 + df_UFDs.log10rho0_u)) * u"Msun/pc^3"
    df_UFDs.rho0_d = (10 .^ (df_UFDs.log10rho0 + df_UFDs.log10rho0_d)) * u"Msun/pc^3"
    return df_UFDs
end