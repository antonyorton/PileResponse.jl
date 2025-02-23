module PileResponse

export
    get_load_vs_depth,
    find_cpt_column_names,
    get_least_squares_interpolator,
    list_available_pile_types,
    get_pile_head_displacement,
    get_initial_pile_head_stiffness,
    get_soil_type_CPT2012,
    get_kc_base_CPT2012,
    get_ultimate_shaft_load,
    get_ultimate_shaft_resistance,
    get_fsmax_shaft_CPT2012,
    get_alpha_shaft_CPT2012,
    get_fsol_shaft_CPT2012,
    get_average_qc_at_pile_base,
    get_E0,
    get_Vs,
    get_Ic,
    get_bigQt,
    get_Fr,
    get_Rf,
    get_qn,
    get_qt,
    get_sigmav0_effective,
    get_sigmav0_total,
    read_delimited_text_file,
    helloworld



include("functions.jl")




function helloworld()
    print("Hi there PileResponse.js, let's see if it works.\n")
    return 123
end

end
