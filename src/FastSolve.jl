"""
module **FastSolve**\\
With this module you can:\\
     - calculate the pile ultimate load in one step via `get_pile_ultimate_load_MN()`\\
     - calculate the pile load response in one step via `get_pile_load_displacement()`
"""
module FastSolve

using Printf
import PileResponse as prs


export
    helloworld,
    get_pile_ultimate_load_MN,
    get_pile_load_displacement


function helloworld(name::AbstractString)
    print("hello there, you must be Mr ", name)
end

"""
    get_pile_load_displacement(cpt_datafile::AbstractString, pile_length::Float64, pile_diameter::Float64; Epile::Int64=30000, pile_type::String="Driven pile - steel closed ended", gamma_soil::Float64=18.5, gw_depth::Float64=3.0, Poisson_ratio::Float64=0.3)

Return `[load (MN), displacement (m)]`.

Allowable pile types:\\
-----------------------------------------\\
 "Driven pile - cast in place"\\
 "Driven pile - steel closed ended"\\
 "Driven pile - pre-cast concrete"\\
 "Driven pile - steel H pile"\\
 "Bored pile - with slurry"\\
 "Bored pile - permanent casing"\\
 "Bored pile - recoverable casing"\\
 "Screw pile - with casing"\\
 "Bored pile - no support"\\
 "Bored pile - with slurry and grooved sockets"\\
 "Driven pile - steel open ended"\\
 "Bored pile - dry bored pile"\\
 "Screw pile - cast in place"\\
 "CFA pile"\\
 "Driven pile - concrete coated steel"\\
"""
function get_pile_load_displacement(cpt_datafile::AbstractString, pile_length::Float64, pile_diameter::Float64; Epile::Int64=30000, pile_type::String="Driven pile - steel closed ended", gamma_soil::Float64=18.5, gw_depth::Float64=3.0, Poisson_ratio::Float64=0.3)

    data = prs.read_delimited_text_file(cpt_datafile, delim=',')

    depth_col, qc_col, fs_col, u2_col =
        prs.find_cpt_column_names(collect(keys(data)))

    depth_m = data[depth_col]
    qc_MPa = data[qc_col]
    fs_MPa = data[fs_col] * 0.001
    u2_MPa = data[u2_col] * 0.001


    # Get Ic, Vs and E0
    # Ic = prs.get_Ic(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73)
    Vs = prs.get_Vs(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73)
    E0 = prs.get_E0(Vs, gamma=gamma_soil, ν=Poisson_ratio)

    # Least squares approx to E0
    fun_E0_linear = prs.get_least_squares_interpolator(depth_m, E0)
    # Get E0 (MPa) at base of pile and at midpoint depth
    E_L = fun_E0_linear(pile_length)
    E_Lon2 = fun_E0_linear(pile_length / 2)

    # Get the initial pile head stiffness k0 (MN/m)
    k0 = prs.get_initial_pile_head_stiffness(pile_length, pile_diameter, Epile, E_L, E_Lon2, ν=Poisson_ratio)

    # # results
    # @printf "E_L = %5.3f MPa\n" E_L
    # @printf "E_Lon2 = %5.3f MPa\n" E_Lon2
    # @printf "k0 = %5.3f MN/m\n" k0

    # Get pile ultimate load
    pile_ult_load_MN = get_pile_ultimate_load_MN(cpt_datafile, pile_length, pile_diameter, Epile=Epile, pile_type=pile_type, gamma_soil=gamma_soil, gw_depth=gw_depth, Poisson_ratio=Poisson_ratio)

    # Pile head loads
    pile_head_loads = 0.01:0.001:0.95*pile_ult_load_MN
    # Pile head displacement
    displacement = prs.get_pile_head_displacement(k0, pile_head_loads, pile_ult_load_MN)

    return stack([pile_head_loads, displacement])
end


"""
    get_pile_ultimate_load_MN(cpt_datafile::AbstractString, pile_length::Float64, pile_diameter::Float64; Epile::Int64=30000, pile_type::String="Driven pile - steel closed ended", gamma_soil::Float64=18.5, gw_depth::Float64=3.0, Poisson_ratio::Float64=0.3)

Return the ultimate pile load in MN.

Allowable pile types:\\
-----------------------------------------\\
 "Driven pile - cast in place"\\
 "Driven pile - steel closed ended"\\
 "Driven pile - pre-cast concrete"\\
 "Driven pile - steel H pile"\\
 "Bored pile - with slurry"\\
 "Bored pile - permanent casing"\\
 "Bored pile - recoverable casing"\\
 "Screw pile - with casing"\\
 "Bored pile - no support"\\
 "Bored pile - with slurry and grooved sockets"\\
 "Driven pile - steel open ended"\\
 "Bored pile - dry bored pile"\\
 "Screw pile - cast in place"\\
 "CFA pile"\\
 "Driven pile - concrete coated steel"\\
"""
function get_pile_ultimate_load_MN(cpt_datafile::AbstractString, pile_length::Float64, pile_diameter::Float64; Epile::Int64=30000, pile_type::String="Driven pile - steel closed ended", gamma_soil::Float64=18.5, gw_depth::Float64=3.0, Poisson_ratio::Float64=0.3)

    data = prs.read_delimited_text_file(cpt_datafile, delim=',')

    depth_col, qc_col, fs_col, u2_col =
        prs.find_cpt_column_names(collect(keys(data)))

    depth_m = data[depth_col]
    qc_MPa = data[qc_col]
    fs_MPa = data[fs_col] * 0.001
    u2_MPa = data[u2_col] * 0.001


    # Get Ic, Vs and E0
    Ic = prs.get_Ic(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73)
    Vs = prs.get_Vs(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73)
    E0 = prs.get_E0(Vs, gamma=gamma_soil, ν=Poisson_ratio)


    # Ultimate shaft load
    # Get soil type
    soil_type_CPT2012 = prs.get_soil_type_CPT2012(Ic)
    prs.get_alpha_shaft_CPT2012()[pile_type]
    # Get fshaft_MPa
    fshaft_MPa = prs.get_ultimate_shaft_resistance(qc_MPa, Ic, pile_type, factor=1.0)
    # ultimate shaft load
    ult_shaft_MN = prs.get_ultimate_shaft_load(depth_m, fshaft_MPa, pile_diameter, pile_length)

    # Ultimate base load
    # Get qc-avg
    qc_avg_base = prs.get_average_qc_at_pile_base(depth_m, qc_MPa, pile_length, pile_diameter, clip_to_30pct=false)
    # Get kc factor
    kc_at_base = prs.get_kc_base_CPT2012()[pile_type][soil_type_CPT2012[depth_m.==pile_length]][1]
    # Calculate ultimate base load
    fb_MPa = kc_at_base * qc_avg_base
    ult_base_MN = pi * pile_diameter^2 / 4 * fb_MPa

    # Pile ulitmate load
    pile_ult_load_MN = ult_base_MN + ult_shaft_MN

    # # Print results
    # @printf "qc_avg_base = %5.3f MPa\n" qc_avg_base
    # @printf "Ultimate shaft load = %5.3f MN\n" ult_shaft_MN
    # @printf "Ultimate base load = %5.3f MN\n" ult_base_MN
    # @printf "Pile ultimate load = %5.3f MN\n" pile_ult_load_MN

    return pile_ult_load_MN
end




end

