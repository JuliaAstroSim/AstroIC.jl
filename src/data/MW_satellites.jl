"""
$(TYPEDSIGNATURES)

Load `src/data/MW_satellites.csv` with units in Galactocentric coordinates.

Reference: Table B.1 and Table B.2 from *battaglia2022gaia_Astron. Astrophys. - Gaia early DR3 systemic motions of local group dwarf galaxies and orbital properties with a massive Large Magellanic Cloud*
"""
function load_data_MW_satellites()
    df_MW_satellites = DataFrame(CSV.File(joinpath(@__DIR__, "MW_satellites.csv"), delim=",", header=true));
    df_MW_satellites.RA = df_MW_satellites.RA * u"°"
    df_MW_satellites.DEC = df_MW_satellites.DEC * u"°"
    df_MW_satellites.Distance = df_MW_satellites.Distance * u"kpc"
    df_MW_satellites.Distance_u = df_MW_satellites.Distance_u * u"kpc"
    df_MW_satellites.Distance_d = df_MW_satellites.Distance_d * u"kpc"
    df_MW_satellites.v_los = df_MW_satellites.v_los * u"km/s"
    df_MW_satellites.v_los_u = df_MW_satellites.v_los_u * u"km/s"
    df_MW_satellites.v_los_d = df_MW_satellites.v_los_d * u"km/s"
    df_MW_satellites.mu_alpha = df_MW_satellites.mu_alpha * u"mas/yr"
    df_MW_satellites.mu_alpha_u = df_MW_satellites.mu_alpha_u * u"mas/yr"
    df_MW_satellites.mu_alpha_d = df_MW_satellites.mu_alpha_d * u"mas/yr"
    df_MW_satellites.mu_delta = df_MW_satellites.mu_delta * u"mas/yr"
    df_MW_satellites.mu_delta_u = df_MW_satellites.mu_delta_u * u"mas/yr"
    df_MW_satellites.mu_delta_d = df_MW_satellites.mu_delta_d * u"mas/yr"

    df_MW_satellites.X = df_MW_satellites.X * u"kpc"
    df_MW_satellites.Y = df_MW_satellites.Y * u"kpc"
    df_MW_satellites.Z = df_MW_satellites.Z * u"kpc"
    df_MW_satellites.v_X = df_MW_satellites.v_X * u"km/s"
    df_MW_satellites.v_Y = df_MW_satellites.v_Y * u"km/s"
    df_MW_satellites.v_Z = df_MW_satellites.v_Z * u"km/s"
    return df_MW_satellites
end