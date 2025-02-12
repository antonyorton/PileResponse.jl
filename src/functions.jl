using DataInterpolations
using DelimitedFiles
using LsqFit


"""
	get_least_squares_interpolator(xvals::AbstractVector, yvals::AbstractVector; p0::AbstractVector = [0.0,0.0])

	Returns a function `fun()` such that `fun(x) = mx + b` is the least squares linear approximation to [xvals, yvals]\n
`p0 [Float64, Float64]` is the initial guess, and `[0.0, 0.0]` seems to work
"""
function get_least_squares_interpolator(xvals::AbstractVector, yvals::AbstractVector; p0::AbstractVector=[0.0, 0.0])
    linearmodel(x, p) = p[1] * x .+ p[2] # define a linear model y = mx + b,  (p = [m,b])
    linearfit = LsqFit.curve_fit(linearmodel, xvals, yvals, p0) #least squares fit to E0
    fun_linear(x) = linearfit.param[1] * x .+ linearfit.param[2] #linear function
    return fun_linear
end


"""
	show_table(columns::AbstractArray, header::AbstractVector{String}; num_rows::Int64 = 5, formatters = ft_printf("%5.3f"))

	Show a formatted table in HTML with given column data and headers. 

`columns` = `[col1data, col2data, ..]` a list of vectors with the data for each column\\
`header` = `["col1 name", "col2 name", ...]` a list of strings with the column names
"""
function show_table(columns::AbstractArray, header::AbstractVector{String}; num_rows::Int64=5, formatters=ft_printf("%5.3f"))
    tabledata = stack(columns)
    n_total = length(tabledata[:, 1])
    indices = Int64.(1.0:floor(n_total / (num_rows - 1)):n_total+1)
    indices[end] = min.(n_total, indices[end])
    return pretty_table(HTML, tabledata[indices, :], formatters=formatters, header=header)
end


"""
	get_pile_head_displacement(k0::Float64, pile_head_loads::AbstractVector{Float64}, pile_ult_resistance::Float64)\n
	Returns the pile head displacement for the given load vector

The assumption is that the pile head stiffness k for a given load P can be approximated as ``k = k_{0}(1 - (P/P_{ult})^{0.3})``.\n
For further details see:\n
 - Mayne, P. W. (2001) Stress-strain-strength flow parameters from enhanced in-situ tests.\n 
 - Fahey, M. and Carter, J. P. (1993) A finite element study of the pressuremeter in sand using a nonlinear elastic plastic model.
\n
`k0` is the initial, or small strain, pile head stiffness\\
`pile_head_loads` is a vector of the load applied at the pile head. It will be sorted into increasing order.\\
`pile_ult_resistance` is the assessed ultimate load capacity of the pile

"""
function get_pile_head_displacement(k0::Float64, pile_head_loads::AbstractVector{Float64}, pile_ult_resistance::Float64)

    # function defining the reduction of k with load
    kdecay(u) = 1 - u^0.3

    pvals = append!([0.0], sort(pile_head_loads))
    pmid = 0.5 * (pvals[2:end] .+ pvals[1:end-1]) #midpoint loads

    kvals = k0 * kdecay.(pmid / pile_ult_resistance)

    displacement = pvals[2:end] ./ kvals

    return displacement
end


"""
	get_initial_pile_head_stiffness(pile_length::Float64, pile_diameter::Float64, Epile::Int64, Esoil_L::Float64, Esoil_Lon2::Float64; ν::Float64 = 0.3)\n

	Returns the initial, or small strain, pile head stiffness

The theory is based on the closed form elastic solution provided by Randolph and Wroth (1978), *Analysis of deformations of vertically loaded piles.*\n
The theory assumes that the soil has a linearly increasing elastic modulus with depth.

`Epile` is the elastic modulus of the pile\\
`Esoil_L` is the small strain (E₀) elastic modulus of the soil at the base of the pile shaft\\
`Esoil_Lon2` is the small strain (E₀) elastic modulus of the soil at the midpoint of the pile shaft\\
`ν` (input as \\nu[tab]) is the Poisson's ratio of the soil

"""
function get_initial_pile_head_stiffness(pile_length::Float64, pile_diameter::Float64, Epile::Int64, Esoil_L::Float64, Esoil_Lon2::Float64; ν::Float64=0.3)

    L = pile_length
    G_L = Esoil_L / 2 / (1 + ν)
    G_Lon2 = Esoil_Lon2 / 2 / (1 + ν)
    r₀ = pile_diameter / 2

    ρ = G_Lon2 / G_L
    ζ = log(2.5 * L * ρ * (1 - ν) / r₀)
    λ = Epile / G_L
    μL = sqrt(2 * L^2 / (ζ * λ * r₀^2))

    numerator = 4 / (1 - ν) + 2 * pi * ρ * L * tanh(μL) / (ζ * r₀ * μL)
    denominator = 1 + 4 / (1 - ν) * L * tanh(μL) / (pi * λ * r₀ * μL)
    k0 = G_L * r₀ * numerator / denominator #Pile head stiffness (MN / m)

    return k0
end


"""
	get_soil_type_CPT2012(Ic::AbstractVector{Float64})\n

	Returns vector of soil types for use in CPT2012 capacity assessment:
		1 = silts and clays (Ic = 2.60 - 3.0)
		2 = intermediate soil (Ic = 2.05 - 2.60)
		3 = sands (Ic = 0.0 - 2.05)

"""
function get_soil_type_CPT2012(Ic::AbstractVector{Float64})

    soil_type = Ic[:]
    soil_type[:] .= 3 # sands (Ic = 0.0 - 2.05)
    soil_type[Ic.>2.05] .= 2 # intermediate soil (Ic = 2.05 - 2.60)
    soil_type[Ic.>2.60] .= 1 # silts and clays (Ic = 2.60 - 3.0)

    return Int.(soil_type)
end


"""
    get_kc_base_CPT2012()\n
	Returns a dictionary with pile base resistance kc factors such that Qult_base = kc * qca, where qca is the equivalent average cone resistance at base (see: get_average_qc_at_pile_base()).\n
`Keys`: pile types\\
`Values`: kc factors for [silt/clay, intermediate soils, sands]

Reference: Frank. R. (2017) Some aspects of pile design in France
"""
function get_kc_base_CPT2012()

    kcvals = Dict()
    kcvals["Bored pile - no support"] = [0.4, 0.3, 0.2]
    kcvals["Bored pile - with slurry"] = [0.4, 0.3, 0.2]
    kcvals["Bored pile - permanent casing"] = [0.4, 0.3, 0.2]
    kcvals["Bored pile - recoverable casing"] = [0.4, 0.3, 0.2]
    kcvals["Bored pile - dry bored pile"] = [0.4, 0.3, 0.2]
    kcvals["Bored pile - with slurry and grooved sockets"] = [0.4, 0.3, 0.2]
    kcvals["CFA pile"] = [0.45, 0.3, 0.25]
    kcvals["Screw pile - cast in place"] = [0.5, 0.5, 0.5]
    kcvals["Screw pile - with casing"] = [0.5, 0.5, 0.5]
    kcvals["Driven pile - pre-cast concrete"] = [0.45, 0.4, 0.4]
    kcvals["Driven pile - concrete coated steel"] = [0.45, 0.4, 0.4]
    kcvals["Driven pile - cast in place"] = [0.45, 0.4, 0.4]
    kcvals["Driven pile - steel closed ended"] = [0.45, 0.4, 0.4]
    kcvals["Driven pile - steel open ended"] = [0.35, 0.3, 0.25]
    kcvals["Driven pile - steel H pile"] = [0.4, 0.4, 0.4]

    return kcvals
end


"""
	get_ultimate_shaft_resistance(qc_MPa::AbstractVector{Float64}, Ic::AbstractVector{Float64}, pile_type::String; factor::Float64 = 1.0)\n

	Returns the shaft resistance (MPa) for each element.

Set `factor` = 1.0 for ultimate resistance, and `factor` < 1.0 for resistances below the ultimate load
"""
function get_ultimate_shaft_resistance(qc_MPa::AbstractVector{Float64}, Ic::AbstractVector{Float64}, pile_type::String; factor::Float64=1.0)
    if (factor > 1.0) || (factor < 0.0)
        return DomainError(factor, "[factor] must be between 0.0 and 1.0 inclusive")
    end

    fsol = get_fsol_shaft_CPT2012(qc_MPa, Ic) .* 0.001
    alpha = get_alpha_shaft_CPT2012()[pile_type]
    fsmax = get_fsmax_shaft_CPT2012()[pile_type] .* 0.001
    soil_type_CPT2012 = get_soil_type_CPT2012(Ic)
    fshaft_unclipped = factor * fsol .* alpha[soil_type_CPT2012]
    # Clip to fsmax
    fshaft = fshaft_unclipped[:]
    fshaft[fshaft_unclipped.>fsmax[soil_type_CPT2012]] = fsmax[soil_type_CPT2012][fshaft_unclipped.>fsmax[soil_type_CPT2012]]

    return fshaft
end


"""
    get_fsmax_shaft_CPT2012()\n
    Returns a dictionary with maximum pile shaft resistance values fsmax,
		such that `Qult_shaft <= fsmax`.\n
`Keys`: pile types\\
`Values`: maximum shaft resistance fsmax(kPa) for [silt/clay, intermediate soils, sands]

Reference: Frank. R. (2017) Some aspects of pile design in France
"""
function get_fsmax_shaft_CPT2012()

    fsmaxvals = Dict()
    fsmaxvals["Bored pile - no support"] = [90.0, 90.0, 90.0]
    fsmaxvals["Bored pile - with slurry"] = [90.0, 90.0, 90.0]
    fsmaxvals["Bored pile - permanent casing"] = [50.0, 50.0, 50.0]
    fsmaxvals["Bored pile - recoverable casing"] = [90.0, 90.0, 90.0]
    fsmaxvals["Bored pile - dry bored pile"] = [90.0, 90.0, 90.0]
    fsmaxvals["Bored pile - with slurry and grooved sockets"] = [90.0, 90.0, 90.0]
    fsmaxvals["CFA pile"] = [90.0, 90.0, 170.0]
    fsmaxvals["Screw pile - cast in place"] = [130.0, 130.0, 200.0]
    fsmaxvals["Screw pile - with casing"] = [50.0, 50.0, 90.0]
    fsmaxvals["Driven pile - pre-cast concrete"] = [130.0, 130.0, 130.0]
    fsmaxvals["Driven pile - concrete coated steel"] = [170.0, 170.0, 260.0]
    fsmaxvals["Driven pile - cast in place"] = [90.0, 90.0, 130.0]
    fsmaxvals["Driven pile - steel closed ended"] = [90.0, 90.0, 90.0]
    fsmaxvals["Driven pile - steel open ended"] = [90.0, 90.0, 50.0]
    fsmaxvals["Driven pile - steel H pile"] = [90.0, 90.0, 130.0]

    return fsmaxvals
end


"""
    get_alpha_shaft_CPT2012()\n
    Returns a dictionary with pile shaft resistance α factors,
		such that Qult_shaft = α * fsol (see: get_fsol_shaft_CPT2012()).\n
`Keys`: pile types\\
`Values`: alpha factors for [silt/clay, intermediate soils, sands]

Reference: Frank. R. (2017) Some aspects of pile design in France
"""
function get_alpha_shaft_CPT2012()

    # Dictionary setting alpha for each pile type, where the entries
    # are for soil types [silt/clay, intermediate soils, sands]

    alphavals = Dict()
    alphavals["Bored pile - no support"] = [0.55, 0.65, 0.7]
    alphavals["Bored pile - with slurry"] = [0.65, 0.8, 1]
    alphavals["Bored pile - permanent casing"] = [0.35, 0.4, 0.4]
    alphavals["Bored pile - recoverable casing"] = [0.65, 0.8, 1]
    alphavals["Bored pile - dry bored pile"] = [0.7, 0.85, 0.9]
    alphavals["Bored pile - with slurry and grooved sockets"] = [0.7, 0.85, 0.9]
    alphavals["CFA pile"] = [0.75, 0.9, 1.25]
    alphavals["Screw pile - cast in place"] = [0.95, 1.15, 1.45]
    alphavals["Screw pile - with casing"] = [0.3, 0.35, 0.4]
    alphavals["Driven pile - pre-cast concrete"] = [0.55, 0.65, 1]
    alphavals["Driven pile - concrete coated steel"] = [1, 1.2, 1.45]
    alphavals["Driven pile - cast in place"] = [0.6, 0.7, 1]
    alphavals["Driven pile - steel closed ended"] = [0.4, 0.5, 0.85]
    alphavals["Driven pile - steel open ended"] = [0.6, 0.7, 0.5]
    alphavals["Driven pile - steel H pile"] = [0.55, 0.65, 0.7]

    return alphavals
end


"""
    get_fsol_shaft_CPT2012(qc_MPa::AbstractVector{Float64}, Ic::AbstractVector{Float64})\n
    Returns fsol(kPa)::Vector the unfactored ultimate pile shaft resistance.

`Ic` is CPT soil behaviour type index

Reference: Frank. R. (2017) Some aspects of pile design in France
"""
function get_fsol_shaft_CPT2012(qc_MPa::AbstractVector{Float64}, Ic::AbstractVector{Float64})

    ysand = 125 .* (1 .- exp.(-0.13 .* qc_MPa)) # Curve for sands, Ic < 2.05
    yclay = 138 .* (1 .- exp.(-0.21 .* qc_MPa)) # Curve for silts and clays Ic > 2.60

    # Clip Ic to between 2.05 and 2.6, with sands below, and silts and clays above
    Ic_clipped = Ic[:] # take a copy of Ic
    Ic_clipped[Ic.<2.05] .= 2.05 # lower (sand) limit
    Ic_clipped[Ic.>2.60] .= 2.60 # upper (silts and clays) limit

    # scaling ratio
    sand_clay_ratio = (Ic_clipped .- 2.05) ./ 0.55  # 0 for sand 1 for clay

    return (1 .- sand_clay_ratio) .* ysand .+ sand_clay_ratio .* yclay
end


"""
	get_average_qc_at_pile_base(depth_m::AbstractVector{Float64},qc_MPa::AbstractVector{Float64}, pile_toe_depth::Float64, pile_diameter::Float64; scale_to_30pct::Bool = false)\n
	Returns average qc within +/- 1.5 pile diameters from the toe, with values scaled to within +/- 30% of the value at the toe

if `clip_to_30pct` = `false`, values will not be scaled to +/- 30% of the value at the toe prior to averaging
"""
function get_average_qc_at_pile_base(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, pile_toe_depth::Float64, pile_diameter::Float64; clip_to_30pct::Bool=false)

    # Get qc values within +/- 1.5 pile diameters from the toe
    base_qcvals = qc_MPa[depth_m.>=pile_toe_depth-1.5*pile_diameter.&&depth_m.<=pile_toe_depth+1.5*pile_diameter]

    # Get qc value at the toe
    base_qc = qc_MPa[argmin(abs.(depth_m .- pile_toe_depth))]


    # limit base_qc values to 0.7 * qc at toe < x < 1.3 * qc at toe
    if clip_to_30pct
        for i in eachindex(base_qcvals)
            if base_qcvals[i] > 1.3 * base_qc
                base_qcvals[i] = 1.3 * base_qc
            elseif base_qcvals[i] < 0.7 * base_qc
                base_qcvals[i] = 0.7 * base_qc
            end
        end
    end

    return mean(base_qcvals)
end


"""
	get_E0(Vs::AbstractVector{Float64}; gamma::Float64 = 18.0, ν::Float64 = 0.3)\n
	returns `E₀(MPa) = 2 * (1 + ν) * gamma / 9.81 * Vs² * 0.001`
`gamma` is soil unit weight in kN/m², assumed constant over `depth_m`\\
`ν` is the Poisson's ratio
"""
function get_E0(Vs::AbstractVector{Float64}; gamma::Float64=18.0, ν::Float64=0.3)
    return 2 * (1 + ν) * gamma / 9.81 * Vs .^ 2 * 0.001
end


"""
	get_Vs(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64 = 18.0, a::Float64 = 0.73)\n
	returns `(αᵥₛ * qn / pa)^0.5`
		where `αᵥₛ = 10^(0.55 * Ic + 1.68)`
`gw_depth` is depth to groundwater in metres\\
`gamma` is soil unit weight in kN/m², assumed constant over `depth_m`\\
`a` is the net area ratio of the cone, typically between 0.70 and 0.85\\
"""
function get_Vs(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0, a::Float64=0.73)
    pa = 0.101 #atmospheric pressure (MPa)
    qn = get_qn(depth_m, qc_MPa, u2_MPa, gamma=gamma, a=a)
    Ic = get_Ic(depth_m, qc_MPa, u2_MPa, gw_depth, gamma=gamma, a=a)
    αvs = 10.0 .^ (0.55 .* Ic .+ 1.68)
    return (αvs .* qn / pa) .^ 0.5
end


"""
	get_Ic(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64 = 18.0, a::Float64 = 0.73)\\
	returns `[(3.47 - log10(Qt))² + (log10(Fr) + 1.22)²]^0.5`
`gw_depth` is depth to groundwater in metres\\
`gamma` is soil unit weight in kN/m², assumed constant over `depth_m`\\
`a` is the net area ratio of the cone, typically between 0.70 and 0.85
"""
function get_Ic(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0, a::Float64=0.73)
    Qt = get_bigQt(depth_m, qc_MPa, u2_MPa, gw_depth, gamma=gamma, a=a)
    Fr = get_Fr(depth_m, qc_MPa, u2_MPa, gamma=gamma, a=a)
    Fr[Fr.<0] .= 0.0
    Qt[Qt.<0] .= 0.0
    Ic = ((3.47 .- log10.(Qt)) .^ 2 .+ (log10.(Fr) .+ 1.22) .^ 2) .^ 0.5
    # Deal with troublesome cases
    Ic[.!isfinite.(Ic)] .= 3.0
    Ic[Ic.>3.0] .= 3.0

    return Ic
end


"""
	get_bigQt(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gw_depth::Float64, gamma::Float64 = 18.0, a::Float64 = 0.73)\\
	returns `Qₜ = (qₜ - σᵥ₀) / σ'ᵥ₀`\\

`gw_depth` is depth to groundwater in metres\\
`gamma` is soil unit weight in kN/m², assumed constant over `depth_m`\\
`a` is the net area ratio of the cone, typically between 0.70 and 0.85
"""
function get_bigQt(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0, a::Float64=0.73)
    qt = get_qt(qc_MPa, u2_MPa, a=a)
    sigmav0 = get_sigmav0_total(depth_m, gamma=gamma)
    sigmav0effective = get_sigmav0_effective(depth_m, gw_depth, gamma=gamma)
    return (qt .- sigmav0) ./ sigmav0effective
end


"""
	function get_Fr(depth_m::AbstractVector{Float64},  qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gamma::Float64 = 18.0, a::Float64 = 0.73)\\
	returns `(fₛ / qₙ) * 100`
`a` is the net area ratio of the cone, typically between 0.70 and 0.85
`gamma` is soil unit weight in kN/m², assumed constant over `depth_m`
"""
function get_Fr(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gamma::Float64=18.0, a::Float64=0.73)
    qn_MPa = get_qn(depth_m, qc_MPa, u2_MPa, gamma=gamma, a=a)
    return fs_MPa ./ qn_MPa * 100
end


"""
	get_Rf(qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; a::Float64 = 0.73)\\
	Returns `(fₛ / qₜ) * 100`

`a` is the net area ratio of the cone, typically between 0.70 and 0.85
"""
function get_Rf(qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; a::Float64=0.73)
    qt_MPa = get_qt(qc_MPa, u2_MPa, a=a)
    return fs_MPa ./ qt_MPa * 100
end


"""
	function get_qn(depth_m::AbstractVector{Float64},  qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gamma::Float64 = 18.0, a::Float64 = 0.73)\\
	returns qₙ = qt - σᵥ₀
`a` is the net area ratio of the cone, typically between 0.70 and 0.85
`gamma` is soil unit weight in kN/m², assumed constant over `depth_m`
"""
function get_qn(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gamma::Float64=18.0, a::Float64=0.73)
    qt_MPa = get_qt(qc_MPa, u2_MPa, a=a)
    sigmav0_total = get_sigmav0_total(depth_m, gamma=gamma)
    return max.(0.0, qt_MPa .- sigmav0_total) #replace negative values with 0.0
end


"""
	get_qt(qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; a::Float64 = 0.73)\\
	Returns `qt = qc + u2(1 - a)`

where `a` is the net area ratio of the cone, typically between 0.70 and 0.85
"""
function get_qt(qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; a::Float64=0.73)
    return qc_MPa .+ u2_MPa .* (1 - a)
end


"""
	get_sigmav0_effective(depth_m::AbstractVector{Float64}, gw_depth::Float64, gamma::Float64 = 18.0)\\
	Returns the effective vertical stress σ'ᵥ₀(MPa) 

`gw_depth` is depth to groundwater in metres\\
`gamma` is soil unit weight in kN/m², assumed constant over `depth_m`
"""
function get_sigmav0_effective(depth_m::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0)

    u_MPa = 0.001 * max.(9.81 * (depth_m .- gw_depth), 0)
    sigma_v0_total = get_sigmav0_total(depth_m, gamma=gamma)

    return sigma_v0_total .- u_MPa
end


"""
	get_sigmav0_total(depth_m::AbstractVector{Float64}; gamma::Float64 = 18.0)\\
	Returns the total vertical stress σᵥ₀(MPa)

`gamma` is unit weight in kN/m²  
"""
function get_sigmav0_total(depth_m::AbstractVector{Float64}; gamma::Float64=18.0)
    return 0.001 * depth_m * gamma
end


"""
    read_delimited_text_file(filepath::String; delim::AbstractChar=',', T::Type=Float64)\\
	Returns a dictionary `data` such that `data[col1_header] = col1_values`

Read a delimited text file. Assumes the first row is the header and that all columns below the header have the same type.
"""
function read_delimited_text_file(filepath::String; delim::AbstractChar=',', T::Type=Float64)
    data, headers = DelimitedFiles.readdlm(filepath, delim, T, header=true)
    return Dict(headers[i] => data[:, i] for i = 1:length(data[1, :]))
end
