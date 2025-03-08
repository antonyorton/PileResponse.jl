using DelimitedFiles
using LsqFit
using Statistics


"""
    get_load_vs_depth(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, Ic::AbstractVector{Float64}, applied_load::Float64, pile_ult_load::Float64, pile_length::Float64, pile_diameter::Float64, pile_type::AbstractString)

Returns `[depth, load]`, which is the `load` carried by the pile at each `depth` for a given `applied_load` at the pile head.

# Example
```julia
mydepth, myload = get_load_vs_depth(depth_m, qc_MPa, Ic, applied_load, pile_ult_load, pile_length, pile_diameter, pile_type)
```
"""
function get_load_vs_depth(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, Ic::AbstractVector{Float64}, applied_load::Float64, pile_ult_load::Float64, pile_length::Float64, pile_diameter::Float64, pile_type::AbstractString)

    load_factor = applied_load / pile_ult_load

    if load_factor > 1.0
        throw(DomainError(applied_load, "applied_load cannot be greater than pile_ult_load"))
    end

    fshaft_factored = get_ultimate_shaft_resistance(qc_MPa, Ic, pile_type, factor=load_factor)

    Qshaft_factored = append!([0.0], pi * pile_diameter * (depth_m[2:end] .- depth_m[1:end-1]) .* 0.5 .* (fshaft_factored[2:end] .+ fshaft_factored[1:end-1]))

    mydepth = depth_m[depth_m.<=pile_length]
    myshaftload = zeros(Float64, length(mydepth))
    myshaftload[1] = load_factor * pile_ult_load
    for i = 2:length(mydepth)
        myshaftload[i] = myshaftload[i-1] - Qshaft_factored[i]
    end

    return [mydepth, myshaftload]
end

"""
    find_cpt_column_names(names::Vector{SubString{String}})

Return a vector `[depth_col, qc_col, fs_col, u2_col]` providing the relevant names for each item from `names`.
    
The function looks for the (non case sensitive) strings "depth", "qc", "fs" and "u2" within the items in `names`.

# Example
```julia
data = read_delimited_text_file("./data/example_cpt_data.csv", delim=',')
column_names = collect(keys(data))
depth_col, qc_col, fs_col, u2_col = find_cpt_column_names(column_names)
```
"""
function find_cpt_column_names(names::Vector{SubString{String}})
    depth_col = [item for item in names if occursin("depth", lowercase(item))][1]
    qc_col = [item for item in names if occursin("qc", lowercase(item))][1]
    fs_col = [item for item in names if occursin("fs", lowercase(item))][1]
    u2_col = [item for item in names if occursin("u2", lowercase(item))][1]
    temp = [depth_col, qc_col, fs_col, u2_col]
    # Return as Vector{String} instead of Vector{SubString{String}}
    return [String(item) for item in temp]
end

"""
	get_least_squares_interpolator(xvals::AbstractVector, yvals::AbstractVector; p0::AbstractVector = [0.0,0.0])

Return a function `fun()` such that `fun(x) = mx + b` is the least squares linear approximation to [xvals, yvals].

- `p0 [Float64, Float64]` is the initial guess, and `[0.0, 0.0]` seems to work.
"""
function get_least_squares_interpolator(xvals::AbstractVector, yvals::AbstractVector; p0::AbstractVector=[0.0, 0.0])
    linearmodel(x, p) = p[1] * x .+ p[2] # define a linear model y = mx + b,  (p = [m,b])
    linearfit = LsqFit.curve_fit(linearmodel, xvals, yvals, p0) #least squares fit to E0
    fun_linear(x) = linearfit.param[1] * x .+ linearfit.param[2] #linear function
    return fun_linear
end


# """
#     show_table(columns::AbstractArray, header::AbstractVector{String}; num_rows::Int64=5, printformat="%5.3f")

# Display a formatted table in HTML with given column data and headers. 

# # Arguments
# - `columns` = `[col1data, col2data, ..]` a list of vectors with the data for each column.
# - `header` = `["col1 name", "col2 name", ...]` a list of strings with the column names.
# """
# function show_table(columns::AbstractArray, header::AbstractVector{String}; num_rows::Int64=5, printformat="%5.3f")
#     tabledata = stack(columns)
#     n_total = length(tabledata[:, 1])
#     indices = Int64.(1.0:floor(n_total / (num_rows - 1)):n_total+1)
#     indices[end] = min.(n_total, indices[end])
#     return pretty_table(HTML, tabledata[indices, :], formatters=ft_printf(printformat), header=header)
# end

"""
    list_available_pile_types()

Return a list of the available pile types which can be used in the analysis.
"""
function list_available_pile_types()
    return sort(collect(keys(get_alpha_shaft_CPT2012())))
end

"""
	get_pile_head_displacement(k0::Float64, pile_head_loads::AbstractVector{Float64}, pile_ult_load::Float64)

Return the pile head displacement for the given load vector.

# Arguments
- `k0::Float64` is the initial pile head stiffness.
- `pile_head_loads::AbstractVector{Float64}` is a vector of the load applied at the pile head. It will be sorted into increasing order.
- `pile_ult_load::Float64` is the assessed ultimate load of the pile.

See also [`get_initial_pile_head_stiffness`](@ref).

The assumption is that the pile head stiffness k for a given load P can be approximated as:

 - ``k = k_{0}(1 - (P/P_{ult})^{0.3})``.


# Reference
 - Mayne, P. W. (2001) Stress-strain-strength flow parameters from enhanced in-situ tests. 
 - Fahey, M. and Carter, J. P. (1993) A finite element study of the pressuremeter in sand using a nonlinear elastic plastic model.

"""
function get_pile_head_displacement(k0::Float64, pile_head_loads::AbstractVector{Float64}, pile_ult_load::Float64)

    # function defining the reduction of k with load
    kdecay(u) = 1 - u^0.3

    pvals = append!([0.0], sort(pile_head_loads))
    pmid = 0.5 * (pvals[2:end] .+ pvals[1:end-1]) #midpoint loads

    kvals = k0 * kdecay.(pmid / pile_ult_load)

    displacement = pvals[2:end] ./ kvals

    return displacement
end


"""
	get_initial_pile_head_stiffness(pile_length::Float64, pile_diameter::Float64, Epile::Int64, Esoil_L::Float64, Esoil_Lon2::Float64; ν::Float64 = 0.3)

Return the initial pile head stiffness for small strain (MN/m).

# Arguments
- `Epile (MPa)` is the elastic modulus of the pile.
- `Esoil_L (MPa)` is the small strain (E₀) elastic modulus of the soil at the base of the pile shaft.
- `Esoil_Lon2 (MPa)` is the small strain (E₀) elastic modulus of the soil at the midpoint of the pile shaft.
- `ν` (input as \\nu[tab]) is the Poisson's ratio of the soil.

The theory is based on the closed form elastic solution provided by Randolph and Wroth (1978). The theory assumes that the soil has a linearly increasing elastic modulus with depth.

# Reference
- Randolph, M. F. and Wroth, C. P. (1978). Analysis of deformation of vertically loaded piles. Jnl. Geot. Eng. Divn., ASCE, 108 (GT12): 1465 - 1488.

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

    # print("rho = ", ρ, "\n")
    # print("zeta = ", ζ, "\n")
    # print("lambda = ", λ, "\n")
    # print("ul = ", μL, "\n")

    numerator = 4 / (1 - ν) + 2 * pi * ρ * L * tanh(μL) / (ζ * r₀ * μL)
    denominator = 1 + 4 / (1 - ν) * L * tanh(μL) / (pi * λ * r₀ * μL)

    # print("numerator = ", numerator, "\n")
    # print("denominator = ", denominator, "\n")
    # print("G_Lr0 = ", G_L * r₀, "\n")

    k0 = G_L * r₀ * numerator / denominator #Pile head stiffness (MN / m)

    return k0
end


"""
	get_soil_type_CPT2012(Ic::AbstractVector{Float64})

Return a vector of soil types (deonted by 1, 2 or 3) for use in CPT2012 capacity assessment:

- 1 = Silts and clays (Ic > 2.70).
- 2 = Intermediate soil (2.50 < Ic <= 2.70).
- 3 = Sands (Ic <= 2.50).
"""
function get_soil_type_CPT2012(Ic::AbstractVector{Float64})

    soil_type = Ic[:]
    soil_type[:] .= 3           # Type 3: sands (Ic <= 2.50)
    soil_type[Ic.>2.50] .= 2    # Type 2: intermediate soil (Ic = 2.50 - 2.70)
    soil_type[Ic.>2.70] .= 1    # Type 1: silts and clays (Ic > 2.70)

    return Int.(soil_type)
end


"""
    get_kc_base_CPT2012()

Return a dictionary with pile base resistance factors ``k_{c}`` such that:

``Q_{ult(base)} = k_{c} \\cdot q_{ca}``, where ``q_{ca}`` is the equivalent average cone resistance at the pile base.

- `Keys`: pile types.
- `Values`: ``k_{c}`` factors for [silt/clay, intermediate soils, sands].

See also [`get_average_qc_at_pile_base`](@ref).

# Reference
- Frank, R. (2017). Some aspects of research and practice for pile design in France. Innov. Infrastruct. Solutions. 2 (32) 1 - 15.
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
    get_ultimate_shaft_load(depth_m::AbstractVector{Float64}, fshaft_MPa::AbstractVector{Float64}, pile_diameter::Float64, pile_length::Float64)

Return the ultimate shaft load for the pile (MN).

`fshaft_MPa` is the ultimate shaft resistance (MPa) for each node corresponding to the `depth_m` vector.
"""
function get_ultimate_shaft_load(depth_m::AbstractVector{Float64}, fshaft_MPa::AbstractVector{Float64}, pile_diameter::Float64, pile_length::Float64)
    ult_shaft_MN = 0.0
    for i in eachindex(depth_m)
        depth_m[i] > pile_length && break
        if i > 1
            ult_shaft_MN += pi * pile_diameter * (depth_m[i] - depth_m[i-1]) *
                            0.5 * (fshaft_MPa[i] + fshaft_MPa[i-1])
        end
    end
    return ult_shaft_MN
end


"""
	get_ultimate_shaft_resistance(qc_MPa::AbstractVector{Float64}, Ic::AbstractVector{Float64}, pile_type::String; factor::Float64 = 1.0)

Return the ultimate shaft resistance (MPa) for each element. The method follows the approach described by Frank (2017) where:
- ``f_{s} = \\alpha\\cdot f_{sol}``, with ``f_{sol} \\leq f_{smax}``
        
Set `factor` = 1.0 for ultimate resistance and `factor` < 1.0 for loads less than the ultimate load.

See also:
- [`get_fsol_shaft_CPT2012`](@ref)
- [`get_alpha_shaft_CPT2012`](@ref)
- [`get_fsmax_shaft_CPT2012`](@ref)

# Reference
- Frank, R. (2017). Some aspects of research and practice for pile design in France. Innov. Infrastruct. Solutions. 2 (32) 1 - 15.
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
    get_fsmax_shaft_CPT2012()
    
Return a dictionary with maximum pile shaft resistance values ``f_{smax}``, such that:

``Q_{ult(shaft)} \\leq f_{smax}``.

- `Keys`: Pile types.
- `Values`: Maximum shaft resistance ``f_{smax}`` (kPa) for [silt/clay, intermediate soils, sands].

# Reference
- Frank, R. (2017). Some aspects of research and practice for pile design in France. Innov. Infrastruct. Solutions. 2 (32) 1 - 15.
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
    get_alpha_shaft_CPT2012()
    
Return a dictionary with pile shaft resistance ``\\alpha`` factors, such that:

``Q_{ult(shaft)} = \\alpha \\cdot f_{sol}``

- `Keys`: Pile types.
- `Values`: ``\\alpha`` factors for [silt/clay, intermediate soils, sands].

See also: [`get_fsol_shaft_CPT2012`](@ref).

# Reference
- Frank, R. (2017). Some aspects of research and practice for pile design in France. Innov. Infrastruct. Solutions. 2 (32) 1 - 15.
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
    get_fsol_shaft_CPT2012(qc_MPa::AbstractVector{Float64}, Ic::AbstractVector{Float64})
    
Return ``f_{sol}`` (kPa), the unfactored ultimate pile shaft resistance.

`Ic` is the CPT soil behaviour type index. If `Ic` is not available, provide a vector of the same length as `qc_MPa` with values 1, 2 or 3 depending on whether the soil is silt/clay, intermediate or sand respectively for each point.

See also [`get_Ic`](@ref)

# Reference
- Frank, R. (2017). Some aspects of research and practice for pile design in France. Innov. Infrastruct. Solutions. 2 (32) 1 - 15.
"""
function get_fsol_shaft_CPT2012(qc_MPa::AbstractVector{Float64}, Ic::AbstractVector{Float64})

    ysand = 125 .* (1 .- exp.(-0.13 .* qc_MPa)) # Curve for sands, Ic < 2.50
    yclay = 138 .* (1 .- exp.(-0.21 .* qc_MPa)) # Curve for silts and clays Ic > 2.70

    # Clip Ic to between 2.50 and 2.70, with sands below, and silts and clays above
    Ic_clipped = Ic[:] # take a copy of Ic
    Ic_clipped[Ic.<2.50] .= 2.50 # lower (sand) limit
    Ic_clipped[Ic.>2.70] .= 2.70 # upper (silts and clays) limit

    # scaling ratio
    sand_clay_ratio = (Ic_clipped .- 2.50) ./ 0.20  # 0 for sand 1 for clay

    return (1 .- sand_clay_ratio) .* ysand .+ sand_clay_ratio .* yclay
end


"""
	get_average_qc_at_pile_base(depth_m::AbstractVector{Float64},qc_MPa::AbstractVector{Float64}, pile_length::Float64, pile_diameter::Float64; clip_to_30pct::Bool = false)

Return the average `qc` within +/- 1.5 pile diameters from the base, with values limited to within +/- 30% of the value at the base.

If `clip_to_30pct` = `false`, values will not be limited to +/- 30% of the value at the pile base prior to averaging.
"""
function get_average_qc_at_pile_base(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, pile_length::Float64, pile_diameter::Float64; clip_to_30pct::Bool=false)

    # Get qc values within +/- 1.5 pile diameters from the toe
    indices = findall(x -> (x >= pile_length - 1.5 * pile_diameter) & (x <= pile_length + 1.5 * pile_diameter), depth_m)
    base_qcvals = qc_MPa[indices]

    # Get qc value at the toe
    base_qc = qc_MPa[argmin(abs.(depth_m .- pile_length))]


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
	get_E0(Vs::AbstractVector{Float64}; gamma::Float64 = 18.0, ν::Float64 = 0.3)
	
Return `E₀(MPa) = 2 * (1 + ν) * gamma / 9.81 * Vs² * 0.001`.

# Arguments
- `Vs` is the shear wave velocity.
- `gamma` is soil unit weight in kN/m², which is assumed constant.
- `ν` is the soil Poisson's ratio, which is assumed constant.

See also [`get_Vs`](@ref).
"""
function get_E0(Vs::AbstractVector{Float64}; gamma::Float64=18.0, ν::Float64=0.3)
    return 2 * (1 + ν) * gamma / 9.81 * Vs .^ 2 * 0.001
end


"""
	get_Vs(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64 = 18.0, a::Float64 = 0.73)

Return `(αᵥₛ * qn / pa)^0.5` where `αᵥₛ = 10^(0.55 * Ic + 1.68)`.

# Arguments
- `gw_depth` is depth to groundwater in metres.
- `gamma` is soil unit weight in kN/m², which is assumed constant.
- `a` is the net area ratio of the cone, typically between 0.70 and 0.85.

# Reference
- Robertson, P. K. and Cabal, K. (2022). Guide to cone penetration testing. 7th Ed. Gregg Drilling LLC.
"""
function get_Vs(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0, a::Float64=0.73)
    pa = 0.101 #atmospheric pressure (MPa)
    qn = get_qn(depth_m, qc_MPa, u2_MPa, gamma=gamma, a=a)
    Ic = get_Ic(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma, a=a)
    αvs = 10.0 .^ (0.55 .* Ic .+ 1.68)
    return (αvs .* qn / pa) .^ 0.5
end


"""
	get_Ic_usingQt_only(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64},
    fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64 = 18.0, a::Float64 = 0.73)

Return the soil behaviour type index based on Qt (not Qtn):
- `Ic = [(3.47 - log10(Qt))² + (log10(Fr) + 1.22)²] ^ 0.5`.

# Arguments
- `gw_depth` is depth to groundwater in metres.
- `gamma` is soil unit weight in kN/m², which is assumed constant.
- `a` is the net area ratio of the cone, typically between 0.70 and 0.85.

# Reference
- Robertson, P. K. and Cabal, K. (2022). Guide to cone penetration testing. 7th Ed. Gregg Drilling LLC.
"""
function get_Ic_usingQt_only(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0, a::Float64=0.73)
    Qt = get_bigQt(depth_m, qc_MPa, u2_MPa, gw_depth, gamma=gamma, a=a)
    Fr = get_Fr(depth_m, qc_MPa, fs_MPa, u2_MPa, gamma=gamma, a=a)
    Fr[Fr.<0] .= 0.0
    Qt[Qt.<0] .= 0.0
    Ic = ((3.47 .- log10.(Qt)) .^ 2 .+ (log10.(Fr) .+ 1.22) .^ 2) .^ 0.5
    # Deal with troublesome cases
    Ic[.!isfinite.(Ic)] .= 3.0
    Ic[Ic.>3.0] .= 3.0

    return Ic
end

"""
	get_Ic(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64},
    fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64 = 18.0, a::Float64 = 0.73)

Return the soil behaviour type index based on Qtn:
- `Ic = [(3.47 - log10(Qtn))² + (log10(Fr) + 1.22)²] ^ 0.5`.

# Arguments
- `gw_depth` is depth to groundwater in metres.
- `gamma` is soil unit weight in kN/m², which is assumed constant.
- `a` is the net area ratio of the cone, typically between 0.70 and 0.85.

# Reference
- Robertson, P. K. and Cabal, K. (2022). Guide to cone penetration testing. 7th Ed. Gregg Drilling LLC.
"""
function get_Ic(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0, a::Float64=0.73)
    Qtn = get_bigQtn(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma, a=a)
    Fr = get_Fr(depth_m, qc_MPa, fs_MPa, u2_MPa, gamma=gamma, a=a)
    Fr[Fr.<0] .= 0.0
    Qtn[Qtn.<0] .= 0.0
    Ic = ((3.47 .- log10.(Qtn)) .^ 2 .+ (log10.(Fr) .+ 1.22) .^ 2) .^ 0.5
    # Deal with troublesome cases
    Ic[.!isfinite.(Ic)] .= 3.0
    Ic[Ic.>3.0] .= 3.0

    return Ic
end


"""
	get_bigQt(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gw_depth::Float64, gamma::Float64 = 18.0, a::Float64 = 0.73)

Return `Qₜ = (qₜ - σᵥ₀) / σ'ᵥ₀`.

# Arguments
- `gw_depth` is depth to groundwater in metres.
- `gamma` is soil unit weight in kN/m², which is assumed constant.
- `a` is the net area ratio of the cone, typically between 0.70 and 0.85.
"""
function get_bigQt(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0, a::Float64=0.73)
    qt = get_qt(qc_MPa, u2_MPa, a=a)
    sigmav0 = get_sigmav0_total(depth_m, gamma=gamma)
    sigmav0effective = get_sigmav0_effective(depth_m, gw_depth, gamma=gamma)
    return (qt .- sigmav0) ./ sigmav0effective
end


"""
	get_bigQtn(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gw_depth::Float64, gamma::Float64 = 18.0, a::Float64 = 0.73)

Return:
- `Qₜₙ = (qₜ - σᵥ₀) / 0.101 * (0.101 / σ'ᵥ₀)ⁿ`.
- where: `n = 0.381 * Ic + 0.05 * σ'ᵥ₀ / 0.101 - 0.15`.

Uses an iterative process.

# Arguments
- `gw_depth` is depth to groundwater in metres.
- `gamma` is soil unit weight in kN/m², which is assumed constant.
- `a` is the net area ratio of the cone, typically between 0.70 and 0.85.
"""
function get_bigQtn(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0, a::Float64=0.73)

    p_atm = 0.101
    qt = get_qt(qc_MPa, u2_MPa, a=a)
    sigmav0 = get_sigmav0_total(depth_m, gamma=gamma)
    sigmav0effective = get_sigmav0_effective(depth_m, gw_depth, gamma=gamma)

    # iteration
    mse_criteria = 0.01 # termination criteria
    mse = 1.0 # mean square error on n
    Fr = get_Fr(depth_m, qc_MPa, fs_MPa, u2_MPa, gamma=gamma, a=a) # get Fr

    # initial guesses
    n = ones(length(depth_m))
    Qtn = (qt .- sigmav0) ./ p_atm .* (p_atm ./ sigmav0effective) .^ n
    Ic = ((3.47 .- log10.(Qtn)) .^ 2 .+ (log10.(Fr) .+ 1.22) .^ 2) .^ 0.5


    for i in (1:40)
        # print("mse = ", round(mse, digits=3), "\n")

        # updated result
        Ic_new = ((3.47 .- log10.(Qtn)) .^ 2 .+ (log10.(Fr) .+ 1.22) .^ 2) .^ 0.5
        n_new = min.(1.0, 0.381 .* Ic_new .+ 0.05 * (sigmav0effective ./ p_atm) .- 0.15)
        Qtn_new = (qt .- sigmav0) ./ p_atm .* (p_atm ./ sigmav0effective) .^ n_new

        # check mse on n
        mse = sum((n_new .- n) .^ 2)
        if mse < mse_criteria
            # print("mse (final) = ", round(mse, digits=3), "\n")
            break
        end

        # save results for next iteration
        Ic = Ic_new
        n = n_new
        Qtn = Qtn_new

    end

    return Qtn
end

"""
	get_Fr(depth_m::AbstractVector{Float64},  qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64},u2_MPa::AbstractVector{Float64}; gamma::Float64 = 18.0, a::Float64 = 0.73)

Return `(fₛ / qₙ) * 100`.

# Arguments
- `gamma` is soil unit weight in kN/m², which is assumed constant.
- `a` is the net area ratio of the cone, typically between 0.70 and 0.85.
"""
function get_Fr(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gamma::Float64=18.0, a::Float64=0.73)
    qn_MPa = get_qn(depth_m, qc_MPa, u2_MPa, gamma=gamma, a=a)
    return fs_MPa ./ qn_MPa * 100
end


"""
	get_Rf(qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; a::Float64 = 0.73)

Return `(fₛ / qₜ) * 100`.

- `a` is the net area ratio of the cone, typically between 0.70 and 0.85.
"""
function get_Rf(qc_MPa::AbstractVector{Float64}, fs_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; a::Float64=0.73)
    qt_MPa = get_qt(qc_MPa, u2_MPa, a=a)
    return fs_MPa ./ qt_MPa * 100
end


"""
	get_qn(depth_m::AbstractVector{Float64},  qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gamma::Float64 = 18.0, a::Float64 = 0.73)

Return `qₙ = qt - σᵥ₀`.

# Arguments
- `gamma` is soil unit weight in kN/m², which is assumed constant.
- `a` is the net area ratio of the cone, typically between 0.70 and 0.85.
"""
function get_qn(depth_m::AbstractVector{Float64}, qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; gamma::Float64=18.0, a::Float64=0.73)
    qt_MPa = get_qt(qc_MPa, u2_MPa, a=a)
    sigmav0_total = get_sigmav0_total(depth_m, gamma=gamma)
    return max.(0.0, qt_MPa .- sigmav0_total) #replace negative values with 0.0
end


"""
	get_qt(qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; a::Float64 = 0.73)
	
Return `qt = qc + u2(1 - a)`.

- `a` is the net area ratio of the cone, typically between 0.70 and 0.85.
"""
function get_qt(qc_MPa::AbstractVector{Float64}, u2_MPa::AbstractVector{Float64}; a::Float64=0.73)
    return qc_MPa .+ u2_MPa .* (1 - a)
end


"""
	get_sigmav0_effective(depth_m::AbstractVector{Float64}, gw_depth::Float64, gamma::Float64 = 18.0)

Return the effective vertical stress `σ'ᵥ₀` (MPa).

# Arguments
- `gw_depth` is depth to groundwater in metres.
- `gamma` is soil unit weight in kN/m², which is assumed constant.
"""
function get_sigmav0_effective(depth_m::AbstractVector{Float64}, gw_depth::Float64; gamma::Float64=18.0)

    u_MPa = 0.001 * max.(9.81 * (depth_m .- gw_depth), 0)
    sigma_v0_total = get_sigmav0_total(depth_m, gamma=gamma)

    return sigma_v0_total .- u_MPa
end


"""
	get_sigmav0_total(depth_m::AbstractVector{Float64}; gamma::Float64 = 18.0)
	
Return the total vertical stress `σᵥ₀` (MPa).

- `gamma` is soil unit weight in kN/m², which is assumed constant.
"""
function get_sigmav0_total(depth_m::AbstractVector{Float64}; gamma::Float64=18.0)
    return 0.001 * depth_m * gamma
end


"""
    read_delimited_text_file(filepath::String; delim::AbstractChar=',', T::Type=Float64)
	
Return a dictionary (`data`) such that:
- `data[col1_header] = col1_values`.

Assumes the first row is the header and that all columns below the header have the same type.

- `delim` is the separator, ',' for comma separated, or '\\t' for tab separated values.
"""
function read_delimited_text_file(filepath::String; delim::AbstractChar=',', T::Type=Float64)
    data, headers = DelimitedFiles.readdlm(filepath, delim, T, header=true)
    return Dict(headers[i] => data[:, i] for i = 1:length(data[1, :]))
end
