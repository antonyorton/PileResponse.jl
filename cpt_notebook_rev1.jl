### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 64afe942-879b-4aea-84ae-cf1fcb766689
# ╠═╡ show_logs = false
begin
    import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.add(url = "https://github.com/antonyorton/PileResponse.jl")
    Pkg.add("PlutoUI")
	Pkg.add("CairoMakie")
	using PileResponse
	using DataInterpolations
	using PlutoUI
	using CairoMakie
	CairoMakie.activate!(type = "svg");
end;


# ╔═╡ ca156a62-ce63-11ef-32b1-bd2ade8dcf4c
begin
	nbsp = @html_str"&nbsp"; #single whitespace for markdown (use $nbsp)	
	nbsp4 = @html_str"&nbsp; &nbsp; &nbsp; &nbsp;"; # 4x whitespace for markdown
end

# ╔═╡ d1b05e7f-a074-43d9-be90-cd8381ea0c39
PlutoUI.TableOfContents(title = "Table of contents")

# ╔═╡ 7427ecd6-0dda-4284-948e-e5902c63daa1
md""" # Pile response based on CPT data
**Author:** Antony Orton $nbsp4 $nbsp4[https://github.com/antonyorton/PileResponse.jl](https://github.com/antonyorton/PileResponse.jl)
"""

# ╔═╡ 38c42433-425f-4f38-a99a-5ed2e4977ac9
md"""Introduction
-----------------
"""

# ╔═╡ 7f303256-ef8d-46e1-b21b-f02d9a4ec0fb
md"""
This notebook is designed to provide an assessment of the load-deflection response of a pile head under static vertical load. The assessment depends on the pile type and the results of a single cone penetrometer test (CPT) conducted to at least the design toe level of the pile.\
\
The assessment of ultimate pile shaft and base capacities follows the approach outlined by Frank (2017). The assessment of the non linear load-deflection response is based on the initial pile head stiffness estimated by the Randolph and Wroth (1978) formula. The pile head stiffness is then factored down for increasing loads following the equation proposed by Mayne (2001).\
\
The interpretation of soil behaviour type (``I_{c}``) and shear wave velocity (``V_{s}``) from the CPT data follow correlations provided by Robertson and Cabal (2022). 

**References**

Frank, R. (2017). Some aspects of research and practice for pile design in France. Innov. Infrastruct. Solutions. 2 (32) 1 - 15.\
\
Fahey, M. and Carter, J. P. (1993). A finite element study of the pressuremeter in sand using a nonlinear elastic plastic model. Can. Geot. Jnl. 30(2): 348 - 362.\
\
Randolph, M. F. and Wroth, C. P. (1978). Analysis of deformation of vertically loaded piles. Jnl. Geot. Eng. Divn., ASCE, 108 (GT12): 1465 - 1488.\
\
Mayne, P. W. (2001). Stress-strain-strength-flow parameters from enhanced in-situ tests. Proc.  Int. Conf. In-Situ Measurement of Soil Properties and Case Histories. Bali, Indonesia, 27 - 48.\
\
Robertson, P. K. and Cabal, K. (2022). Guide to cone penetration testing. 7th Ed. Gregg Drilling LLC
\
\
"""

# ╔═╡ 387f30fb-3b91-405d-bb11-436f3e45f188
md"""Pile test details
--------
A pile load test will be carried out at Mobile, Alabama. The test comprises a 457 mm diameter, 9.5 mm wall, closed-toe, steel-pipe pile to be driven to 17.0 m below ground and grouted after driving. The pile will stick up approximately 0.9 m above ground.\
\
**Site conditions**\
\
The site is level, with the soil profile consisting of material with 90% sand, apart from a zone between 4.5 m and 6.5 m of gravelly sand with 30% fines. The depth to groundwater is 3.0 m.
\
\
"""

# ╔═╡ ab327aef-981f-4898-be56-bfac4852ebd8
md""" Soil and pile properties
--------
"""

# ╔═╡ 7b404570-cebd-44af-8b87-7179cc2fb1f3
md"""
**Soil properties** (assumed constant throughout the profile).
"""

# ╔═╡ d632faf0-75ee-489b-81b4-9b8575a083e2
md"Soil unit weight (kN/m³)"

# ╔═╡ 1e50db38-e0ad-4f88-9baa-8d708d577e5f
gamma_soil = 20.5;

# ╔═╡ ebe80e52-af49-4475-a32f-d61b36ec8716
md"Soil poisson ratio"

# ╔═╡ bf2c5ffc-d7d1-424d-9a4e-0f78b0ae5ebd
Poisson_ratio = 0.3;

# ╔═╡ 8a3d5f67-4e4f-4c33-904c-5bde4a9cc3ef
md"**Pile type**"

# ╔═╡ 2a956b1c-be86-4218-a8d4-e6baafe048dd
@bind pile_type Select(sort(collect(keys(get_alpha_shaft_CPT2012()))))

# ╔═╡ 1b9997e7-a282-453d-b5a7-a4705bc4ad84
md"""
**Allowable pile head settlement (m)**\
(The pile capacity will be the load at which the allowable settlement is reached)
"""

# ╔═╡ 8f3d0cc8-6466-43dd-82c8-4d09348e7fc5
allowable_pile_head_settlement_m = 0.060;

# ╔═╡ ac1b64ea-b0cb-4d81-8aa9-214c3bd3d42f
md"""
**Pile geometry and Young's modulus**
"""

# ╔═╡ b8301989-928b-409f-95bb-ca15fb5a2431
md"Pile toe depth (meters below ground)"

# ╔═╡ db43b118-16b9-4d66-be86-a129fc286107
pile_toe_depth = 17.0;

# ╔═╡ 8435f38b-91a6-47e1-b596-2296e16070f0
md"Pile diameter (metres)"

# ╔═╡ 88c706e2-7e38-4f3f-a257-16d522f964e7
pile_diameter = 0.457;

# ╔═╡ 13a9ad72-1687-4b6c-998e-035be8e93918
md"""Estimation of ``E_{pile}``"""

# ╔═╡ eb7b6073-d5b3-4bde-b390-5cc76b7d2cff
begin
	# Calculation of Epile assuming a steel shell, grout filled pile
	odiam = 0.457
	idiam = 0.457 - 0.0095
	Aouter = pi * odiam ^ 2 / 4
	Ainner = pi * idiam ^ 2 / 4
	Eouter = 200 # Gpa
	Einner = 20  # GPa
	Epile = ((Aouter - Ainner) * Eouter + Ainner * Einner )/ Aouter
end

# ╔═╡ 7e9b0f08-1304-47f6-b981-9a4c001df324
md"Pile Young's modulus (MPa)"

# ╔═╡ 2df2143c-0e08-4bb6-b056-d56cdd4eb21a
Epile_MPa = 27400;

# ╔═╡ 6edf2a4b-4672-437e-b8da-d0d91c9f924f
md"**Groundwater level**"

# ╔═╡ a184c4eb-cb0e-4a99-8568-9b796ad0ab57
md"Depth to groundwater (metres below ground)"

# ╔═╡ 2c47f0aa-22ea-4723-8c85-cadba3443472
gw_depth = 3.0;

# ╔═╡ 805713cf-e4fa-44cd-b694-5ebe05d7637f
md"""
Read CPT data file and assign keys to variables
-------------
In this section we read the CPT data file. The file must be a delimited text file such as a .csv.\
\
Once the file is read into the notebook, it is important that the variables `qc_col`, `depth_col`, `fs_col` and `u2_col` are assigned to the relevant column names (denoted as 'keys') from the input file.\
\
The CPT data must include ``q_{c}``, ``f_{s}`` and ``u_{2}`` data. Note, however, that a similar assessment could be completed with just ``q_{c}`` data and knowledge of the soil type. See Frank (2017) for further information.\
\

"""

# ╔═╡ acb1bb78-6ff9-4b9a-85bd-b56c726635a1
data = read_delimited_text_file("./data/example_cpt_data.csv",delim = ',');

# ╔═╡ d8728f41-547d-4e4e-a1d0-689a53ac7d28
md"""
!!! attention
    Carefully check the (relative) path and the delimeter in the above function.
"""

# ╔═╡ 46790341-9cae-455b-b0b1-231c10f5c6e1
md"Obtain the column names from the input file."

# ╔═╡ a8117f4a-a2fb-4a5f-bce2-0f88819a22a8
column_names = collect(keys(data));

# ╔═╡ 7e38b7c1-5325-4c94-96c8-b06a2fc2c561
md"""**Assign keys to variables** - this is automated, however the input file must have column names which contain "qc", "fs", "u2" and "depth"."""

# ╔═╡ 715576a6-6f90-43c6-9e81-6af467429ae6
begin
	depth_col = [item for item in column_names if occursin("depth",lowercase(item))][1]
	qc_col = [item for item in column_names if occursin("qc",lowercase(item))][1]
	fs_col = [item for item in column_names if occursin("fs",lowercase(item))][1]
	u2_col = [item for item in column_names if occursin("u2",lowercase(item))][1]
end;

# ╔═╡ ec1769aa-2424-4aed-8696-e157fab03673
md"""Review that the following assignments match the input file column names. If not, update your input file to use common names "qc", "fs", "u2" and "depth"."""

# ╔═╡ 779a21a3-82b0-4d22-95dd-d1f896b3a272
[qc_col, depth_col, fs_col, u2_col]

# ╔═╡ a3337224-351b-4776-9226-36ece5318f56
md"""
!!! attention
    Check that the above list corresponds to the correct column names from the input file.
"""

# ╔═╡ ece64ce0-f597-4846-87d5-4146b7277260
md"""Smooth the CPT data
-----------------------------------
"""

# ╔═╡ 811a0b94-aaf9-4467-89e1-2723651dfa1c
md"""
In this section, the user specifies parameters to assign the number of nodes along the pile, and to smooth the CPT data for analysis. The parameters required are:
 - `num_steps_along_pile`: The number of desired nodes along the pile.
 - `knots`: A smoothing parameter where a lower value results in more smoothing.
"""

# ╔═╡ 29115972-de54-48b7-b91e-3b747f249f4a
md"""**Smooth the data** using a B-Spline approximation\
 $(nbsp4) (see [DataInterpolations.jl](https://docs.sciml.ai/DataInterpolations/stable/) for further details)"""

# ╔═╡ e159c66d-36c2-409a-b287-f07e381f8c6a
md"**Parameters** to adjust the smoothing and coarseness of the data for analysis"

# ╔═╡ 171bdb61-6e79-408a-8835-1706a415e688
knots = 45; num_steps_along_pile = 120;

# ╔═╡ 679e440d-2e39-4cb7-a1c9-9c9ceeac9f6d
begin
	fun_qc = DataInterpolations.BSplineApprox(data[qc_col], data[depth_col], 3, knots, :Uniform, :Uniform, extrapolation=ExtrapolationType.Linear);
	fun_fs = DataInterpolations.BSplineApprox(data[fs_col], data[depth_col], 3, knots, :Uniform, :Uniform, extrapolation=ExtrapolationType.Linear);
	fun_u2 = DataInterpolations.BSplineApprox(data[u2_col], data[depth_col], 3, knots, :Uniform, :Uniform, extrapolation=ExtrapolationType.Linear);
end;

# ╔═╡ d6178e74-b9cb-4587-bc58-6bc1f42f9e2e
md"""
!!! warning
    Some values for `knots` may result in error and it is recommended to try a few different values.
"""

# ╔═╡ 12a5793a-30de-4f49-a95a-81bb0e4f48a3
md"**Create a new depth grid** based on the input `num_steps_along_pile`."

# ╔═╡ 3d0b9618-f5ba-4b2e-a559-f5862d05b485
new_depth = collect(0.01:(pile_toe_depth - 0.01)/num_steps_along_pile:data[depth_col][end]);

# ╔═╡ d1024dbc-436e-44ed-96cc-1e030547575c
md""" Use the new grid and smoothed values for the interpretation."""

# ╔═╡ 5f5e2b35-629c-4a11-89a9-4b129894053b
begin
	depth_m = new_depth;
	qc_MPa = fun_qc(new_depth);
	fs_MPa = 0.001 * fun_fs(new_depth);
	u2_MPa = 0.001 * fun_u2(new_depth);
end;

# ╔═╡ c166a7d0-39b8-4f03-bab7-984722cdfda2
md"The adopted `num_steps_along_pile` results in a spacing of $(round(depth_m[2] - depth_m[1], digits=2)) m between nodes"

# ╔═╡ 2c996085-e5b6-4d63-99db-6a47af2cd13b
md"""Obtain ``I_{c}`` and ``V_{s}`` from CPT correlations
-----------------------------------
"""

# ╔═╡ 3968f88c-d0ea-4b6e-aea2-407584ad917e
md"""For details of the correlations below please refer to Robertson and Cabal (2022)."""

# ╔═╡ f7bb1887-3a38-4752-9ee2-769734bed772
md"""
#### Get Ic
The soil behaviour type is given by:

- ``I_{c} = \sqrt{(3.47 - log_{10}(Q_{t}))^{2} + (log_{10}(F_{r}) + 1.22)^{2}}``
"""

# ╔═╡ 8414a932-8391-4e47-9912-f7039f1ccdb3
Ic = get_Ic(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma = gamma_soil, a = 0.73);

# ╔═╡ 7e99bf7f-4b8d-41a2-9319-887244f75a9c
md"""
#### Get Vs
The shear wave velocity is given by:\

- ``\sqrt{(\alpha_{vs}\cdot \large\frac{q_{n}}{\small{0.101}})}`` $(nbsp) where ``\alpha_{vs} = 10^{(0.55 Ic + 1.68)}``
"""

# ╔═╡ 8f23be3d-ba6a-4340-8169-7881eaf7ddea
	Vs = get_Vs(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma = gamma_soil, a = 0.73);

# ╔═╡ 93606b0b-7294-424a-8edc-3e7901fae6fd
md"""
#### Get E₀
The small strain elastic modulus is given by:\

- ``E_{0}\small{(MPa)} = 2 (1 + ν)\cdot \frac{\gamma}{9.81}\cdot V_{s}^{2}\cdot 0.001``, $(nbsp) where ``\gamma`` is the unit weight of soil in kN/m³
"""

# ╔═╡ 075ec3aa-666a-43d5-a551-10e6b9a10d3a
E0 = get_E0(Vs, gamma = gamma_soil, ν = Poisson_ratio);

# ╔═╡ f5546bf5-78e7-413e-b9f3-9bb6005438b6
md"""**Create a linear least squares** fit to ``E_{0}``."""

# ╔═╡ 81a00b17-402d-4a20-a3ab-38b6eef32f6f
fun_E0_linear = get_least_squares_interpolator(depth_m, E0);

# ╔═╡ f32428f3-b5c7-430f-bb2b-f9be098fcf53
md"""Plot CPT data
-----------------
"""


# ╔═╡ c10cfbb3-1f57-4884-a2bc-f2e5827e0b2e
md"""
The smoothed data used in the assessment is shown as a black line and the raw data as green dots.
"""

# ╔═╡ 3de17617-19d6-4550-8b4b-7bd2d38e3aee
md"""Plot assessed ``I_{c}``, ``V_{s}`` and ``E_{0}`` 
--------------------------------
"""

# ╔═╡ a1795fc2-61b1-4ccc-a39a-f3a18a5b38dd
md"""
 - In the ``I_{c}`` plot, the soil types used in the pile assessment are shown in blue, where 1, 2 and 3 represent silts-clays, intermediate soils and sands, respectively.\
 - In the ``E_{0}`` chart, a linear least squares fit is shown as a dashed line.
"""

# ╔═╡ 97ed1519-b12e-4d72-980b-b5fca12d5cdb
md"""Pile ultimate resistance
-----------------------------
"""

# ╔═╡ 5dca4557-ef1f-4737-a66a-c080decb6460
md"""--------------"""

# ╔═╡ e5a03703-8dfb-413a-8df9-67896cf9b2f4
md"""
**Soil types**\
The calculation of pile ultimate resistance, following the appoach of Frank (2017), relies on several factors which require a soil type.
Using the assessed ``I_{c}`` value, the soil types required for the assessment are defined in this notebook as follows:\
\
$nbsp4 ``I_{c} > 2.60`` : Silt and clay (soil type 1).\
$nbsp4 ``2.05 < I_{c} \leq 2.60`` : Intermediate soil (soil type 2).\
$nbsp4 ``I_{c} \leq 2.05`` : Sand (soil type 3).
\
\
The assessed soil types are shown on the chart for ``I_{c}`` above.
"""

# ╔═╡ 00a15171-61c0-4ed8-8db6-203cc288fba6
soil_type_CPT2012 = get_soil_type_CPT2012(Ic);

# ╔═╡ 3cb2c6ba-a463-4e58-9d3b-3d8577363784
md"""
### Ultimate shaft resistance
"""

# ╔═╡ d980942c-e160-472b-98eb-a5a939e04d1a
md"""This calculation is provided in the `get_ultimate_shaft_resistance` function."""

# ╔═╡ d8c6bcc9-355d-46ca-b13a-e130e8895ede
begin
	fshaft_MPa = get_ultimate_shaft_resistance(qc_MPa, Ic, pile_type, factor = 1.0);
	ult_shaft_MN = get_ultimate_shaft_load(depth_m, fshaft_MPa, pile_diameter, pile_toe_depth);
end;

# ╔═╡ 433be51b-8738-42d2-92c0-0f62344aab3f
md"""
The ultimate shaft resistance, **``f_{s}`` = $(round(ult_shaft_MN, digits = 2)) MN**, is obtained as follows:\
\
 $nbsp $nbsp``f_{s} = \alpha\cdot f_{sol}`` $nbsp $nbsp $nbsp $nbsp with ``f_{sol} \leq f_{smax}`` $nbsp where:\

- ``f_{sol}`` is the unfactored shaft resistance dependent on ``q_{c}`` and soil type.
- ``\alpha`` is a factor dependent on soil type and pile type.\
- ``f_{smax}`` is a the limiting shaft resistance dependent on soil type and pile type.
"""

# ╔═╡ fd8d59af-d979-4bcf-8aef-7ed338f95e71
md"""
The ``f_{sol}`` function is shown below for sands and clays. Intermediate soils are somewhere in between.
"""

# ╔═╡ 48be6b46-9061-483e-a574-4870d9f33669
begin
	#Data
	qcfsol = (0:0.1:20);
	fsol_sands = get_fsol_shaft_CPT2012(qcfsol,ones(length(qcfsol)))
	fsol_clays = get_fsol_shaft_CPT2012(qcfsol,3 * ones(length(qcfsol)))
	#Figure
	fig_fsol = Figure(size = (700,400));
	ax_fsol = Axis(fig_fsol[1,1], ylabel = "fsol (kPa)", xlabel = "qc (MPa)", yticks = (0:25:150))
	lsand1 = lines!(ax_fsol, qcfsol, fsol_sands, color = :Red)
	lclay1 = lines!(ax_fsol, qcfsol, fsol_clays, color = :Blue)
	#Legend
	axislegend(ax_fsol,[lsand1, lclay1], ["Sand and gravel", "Clay and silt"], position = :lt);
end;

# ╔═╡ 392f8de0-ae6b-40fe-9b77-cd4603831eec
	fig_fsol

# ╔═╡ 8cd9047c-3796-46f0-a27f-d98ef86888a6
md"""
### Ultimate base resistance
"""

# ╔═╡ 234d9744-0550-4485-82d1-b6a74ba6f1be
qc_avg_base = round(get_average_qc_at_pile_base(depth_m, qc_MPa, pile_toe_depth, pile_diameter, clip_to_30pct = false), digits = 3)

# ╔═╡ 673a87a4-569d-4c33-ae3e-98debc89df73
kc_at_base = get_kc_base_CPT2012()[pile_type][soil_type_CPT2012[depth_m .== pile_toe_depth]][1]

# ╔═╡ d50e9d4c-2962-4735-bea8-4bd20d120f32
begin
	fb_MPa = kc_at_base * qc_avg_base;
	ult_base_MN = pi * pile_diameter ^2 / 4 * fb_MPa
end;

# ╔═╡ a30fb07b-ec4a-4e41-be61-70718894516a
md"""
The ultimate base load of **$(round(ult_base_MN, digits = 2)) MN** is obtained as follows:\

 
 $nbsp $nbsp``f_{b} = k_{c}\cdot q_{ca}`` $nbsp where:\


 $nbsp $nbsp $nbsp $nbsp``k_{c}`` is a factor dependent on soil type and pile class, and;\
$nbsp $nbsp $nbsp $nbsp``q_{ca}`` is the equavalent average cone resistance for the soil within 1.5 pile diameters of the base.
\
\
"""

# ╔═╡ 35138641-da56-41d2-a266-3e0de45a832d
md"The ultimate pile load is then the sum of the ultimate shaft and base resistances."

# ╔═╡ f40fc671-1cd6-4b31-a1b3-a5da433f8805
pile_ult_load = round(ult_shaft_MN + ult_base_MN, digits = 3);

# ╔═╡ b4291248-c6f6-49ae-b702-b64169a96b3a
md"""
The ulitimate resistance of the pile is **$(round(pile_ult_load, digits = 2)) MN**, comprising:

 * Ultimate shaft resistance of $(round(ult_shaft_MN, digits = 2)) MN.
 * Ultimate base resistance of $(round(ult_base_MN, digits = 2)) MN.
The sections below provide details on the caculation.
"""

# ╔═╡ 03e1dd77-669d-46b2-bff0-8e23fef7cf39
md"""Pile design capacity
-------------------------
"""

# ╔═╡ 3a4fc7e3-cf0c-4ada-ae84-17de8ac74f52
md"""
!!! note

	The `allowable_pile_head_settlement_m` can be adjusted in the pile properties section above.
"""

# ╔═╡ b518a716-062b-4e91-bff8-cbd51e09be1d
md"""-----------------"""

# ╔═╡ 15715b84-c935-4903-ae55-dc115c530dab
md"""Pile load-displacement response
---------------------------------
"""

# ╔═╡ 0b104fff-80fe-43e9-be05-99521aeb0535
begin
	E_L = fun_E0_linear(pile_toe_depth)
	E_Lon2 = fun_E0_linear(pile_toe_depth / 2)
	k0 = get_initial_pile_head_stiffness(pile_toe_depth, pile_diameter, Epile_MPa, E_L, E_Lon2, ν = Poisson_ratio)
end;

# ╔═╡ 7d5bfca5-5ebd-4e5a-8d7f-0e9c1684a922
md"""
The small strain elastic modulus along the pile shaft is assumed to vary linearly, and a least squares approximation gives:\
\
 $(nbsp4)``E_{L}``: ( $(floor(Int64, round(E_L))) MPa ) The small strain elastic modulus at the base of the shaft.\
$(nbsp4)``E_{L/2}``: ( $(floor(Int64, round(E_Lon2))) MPa ) The small strain elastic modulus at the midpoint of the shaft.\
"""

# ╔═╡ 64013deb-af90-465d-a8f5-581e0f6baab6
md"""The initial pile head stiffness, ``k_{0}``, taking account of pile compressibility, is computed following the closed form elastic solution by Randolph and Wroth (1978) as:\
\
 $(nbsp4)``k_{0}``: ( $(floor(Int64, round(k0))) MN/m ) The initial pile head stiffness.\
"""

# ╔═╡ fe095cea-69a2-46a1-823c-dc1bb7613469
md"""
### Pile load displacement curve
"""

# ╔═╡ efe4f0cc-b42b-4898-9103-0ca91747d1c8
md"""
The load displacement curve is derived following a method proposed by Mayne (2001), based on work by Fahey and Carter (1993), which assumes that the pile head stiffness varies as a function of the load ratio ``P/P_{ult}``:\
\
 $(nbsp4)``k = k_{0} \cdot (1 - (P/P_{ult})^{0.3})``
\
\
"""

# ╔═╡ e5af32ab-6eb3-46d5-8a13-6d6efb8c8e9d
pile_head_loads = 0.01:0.001:0.90*pile_ult_load; #use small load step

# ╔═╡ 4f9e4011-e02e-4fb7-a76f-d4ac0ca3ad43
displacement = get_pile_head_displacement(k0, pile_head_loads, pile_ult_load);

# ╔═╡ dc029c10-1afb-4401-85c9-732aa7962510
pile_capacity_MN = round(pile_head_loads[argmin(abs.(displacement .- allowable_pile_head_settlement_m))],digits = 2);

# ╔═╡ 380622fb-3e25-4259-9ae0-45af9fd88d0c
md"""
The pile design capacity is **$(pile_capacity_MN) MN**.\
\
This is the load at which the (user defined) allowable pile head settlement of $(allowable_pile_head_settlement_m) m is reached.
"""

# ╔═╡ c232752c-4f7e-4b32-92f0-e793f4838975
displacement[argmin(abs.(displacement .- allowable_pile_head_settlement_m))]

# ╔═╡ 741eb051-a7cc-456a-b992-044d0acd460a
md""" 
**Pile type:** $(pile_type).\
**Pile diameter:** $(pile_diameter) m.\
**Pile toe depth:** $(pile_toe_depth) m.\
"""

# ╔═╡ 4e9ccc3b-e977-4f15-b0c5-7c5daca84901
begin
	indices = pile_head_loads .< pile_capacity_MN
	show_table(
		[pile_head_loads[indices], displacement[indices]],
		["Load (MN)", "Displacement (m)"],
		num_rows = 15,
		printformat = "%8.6f"
	)
end


# ╔═╡ c46f48bb-e03b-4b90-95cb-513039d5416f
begin
	plot_indicies = pile_head_loads .< pile_capacity_MN
	figDisp = Figure(size = (600,400))
	Axis(figDisp[1,1], xticks = (0:0.01:0.2), yticks = (0:0.5:pile_capacity_MN), xlabel = "Displacement (m)", ylabel = "Pile head load (MN)")
	lines!(figDisp[1,1], displacement[plot_indicies], pile_head_loads[plot_indicies])
end;

# ╔═╡ 8f9a73f7-d3ab-4f05-b5cc-7daeab7afd66
figDisp

# ╔═╡ fef5eba2-7457-4bdf-a889-3650ba38335f
md"""### Load carried by pile shaft at capacity
"""

# ╔═╡ 723a758c-941e-450f-ac2c-60b6a7e140db
capacity_factor = round(pile_capacity_MN / pile_ult_load, digits = 3);

# ╔═╡ c58cb696-e3f6-441d-81f1-5e4915db081d
md""" The pile capacity is assessed to be $(round(100 * capacity_factor, digits = 0)) % of the ultimate load."""

# ╔═╡ 9ea12381-d2e4-442a-bb67-95fdf5f67437
md"""Calculate the resistance (MPa) for each element of the shaft"""

# ╔═╡ 97e12fce-2e96-49fd-8476-0ea178de8c97
fshaft_factored = get_ultimate_shaft_resistance(qc_MPa, Ic, pile_type, factor = capacity_factor);

# ╔═╡ 8335941c-dc8a-4dd9-a558-5314c722c523
md"""Calculate the load (MN) carried by each element of the shaft"""

# ╔═╡ 338baa40-0742-4e48-9aca-6e7489c45e1b
Qshaft_factored = append!([0.0],pi * pile_diameter * (depth_m[2:end] .- depth_m[1:end-1]) .* 0.5 .* (fshaft_factored[2:end] .+ fshaft_factored[1:end-1]));

# ╔═╡ 53f5866b-eecc-423a-9f64-b22a0bfc2de9
md"""Calculate the load (MN) carried by the shaft below a given depth"""

# ╔═╡ bb86eff9-f3ea-4ade-ae3c-5feeac9d8e34
begin 
	mydepth = depth_m[depth_m .<= pile_toe_depth]
	myshaftload = zeros(Float64,length(mydepth))
	myshaftload[1] = capacity_factor * pile_ult_load
	for i=2:length(mydepth)
		myshaftload[i] = myshaftload[i-1] - Qshaft_factored[i]
	end
	# get the loads carried by shaft and base at capacity
	total_shaft_MN_at_capacity = pile_ult_load * capacity_factor - myshaftload[end];
	total_base_MN_at_capacity = pile_ult_load * capacity_factor - total_shaft_MN_at_capacity;
 end;

# ╔═╡ 6c5bdbe9-3253-4f6e-87aa-d62ba21016d2
md"""### Load distribution along pile shaft at capacity"""

# ╔═╡ f2a332dd-c5fa-4a89-a2d7-89818f3f23b3
md"""
The load carried by the pile at capacity is $(pile_capacity_MN) MN, comprising:\

 -  $(round(total_shaft_MN_at_capacity, digits = 3)) MN carried by the shaft\
 -  $(round(total_base_MN_at_capacity, digits = 3)) MN carried by the base
"""

# ╔═╡ 789627e2-0e67-4692-b750-a155b9f16790
show_table([-mydepth, myshaftload],["Elevation (m)", "Load (MN)"], num_rows=15, printformat = "%5.3f")

# ╔═╡ 4e9ef2a1-e2a0-4209-8eb5-6212fb4fdb5f
md"""---------------------
\
End of user input
\
"""

# ╔═╡ c0be21b4-a70a-4dc5-b2b9-dba8ff243a99
md"""-------------
## Appendix 1: Plot configurations
"""

# ╔═╡ dc5cb03e-2fda-4f0c-a580-0deada74707d
md"""Common plot parameters"""

# ╔═╡ 650d88ea-5bdd-4fb1-b734-9184774681aa
begin
	ytick_min = -floor(round(depth_m[end], digits = 2)) - 1.0; #ymin
	mylimits = (0.0, nothing, ytick_min, 0.001) #(xmin,xmax,ymin,ymax)
	markersize = 4.0;
	rasterize = 3;
	axislabelsize = 10
	defaultlinecolour = color = RGBAf(0.2,0.6,0.2)
end;

# ╔═╡ b1b1c43f-cc35-4c18-a4ff-5ea0fb8996e6
begin
	figShaftLoad = Figure(size = (350,600))
	Axis(figShaftLoad[1,1], xticks = (0:0.5:pile_ult_load), yticks = (ytick_min:0), limits = ((0,pile_ult_load), (-pile_toe_depth - 1, 0)), xlabel = "Load (MN)", ylabel = "Elevation (m)")
	lines!(myshaftload, -mydepth)
	figShaftLoad
end

# ╔═╡ 2fefead9-8aab-4044-92bb-2892b15836ca
md"Plot configuration for qc, fs and u2"

# ╔═╡ ff57b7dd-a174-4f15-9cd9-1bb08087c083
begin
	#Plot configuration for qc, fs and u2
	fig_cpt = Figure(size = (600,600));
	#Axes
	p1 = fig_cpt[1,1]; ax1 = Axis(p1, ylabel = "Elevation (m)", xlabel = "qc (MPa)", yticks = (ytick_min:0), yticklabelsize = axislabelsize, xticklabelsize = axislabelsize, limits = mylimits);
	p2 = fig_cpt[1,2]; ax2 = Axis(p2, xlabel = "fs (kPa)", yticks = (ytick_min:0.0), yticklabelsize = axislabelsize, xticklabelsize = axislabelsize, limits = mylimits);
	p3 = fig_cpt[1,3]; ax3 = Axis(p3, xlabel = "u2 (kPa)", yticks = (ytick_min:0), yticklabelsize = axislabelsize, xticklabelsize = axislabelsize, limits = (nothing, nothing, ytick_min, 0.001));
	#plots
	scatter!(fig_cpt[1,1], data[qc_col], -data[depth_col], rasterize = rasterize, markersize = markersize, color = defaultlinecolour);
	scatter!(fig_cpt[1,2], data[fs_col], -data[depth_col], rasterize = rasterize, markersize = markersize, color = defaultlinecolour );
	scatter!(fig_cpt[1,3], data[u2_col], -data[depth_col], rasterize = rasterize, markersize = markersize, color = defaultlinecolour );
	#smoothed data
	lines!(fig_cpt[1,1], fun_qc(new_depth), -new_depth, rasterize = rasterize, color=RGBf(0.1, 0.1, 0.1));
	lines!(fig_cpt[1,2], fun_fs(new_depth), -new_depth, rasterize = rasterize, color=RGBf(0.1, 0.1, 0.1));
	lines!(fig_cpt[1,3], fun_u2(new_depth), -new_depth, rasterize = rasterize, color=RGBf(0.1, 0.1, 0.1));
	# save("figure.pdf", fig_cpt, pdf_version="1.4");
end;
	

# ╔═╡ 201ebe2b-fdaa-430d-b1d4-090f4c39a9b7
	fig_cpt

# ╔═╡ 53ef7ea8-93f5-466d-a7ae-66819ac8c1dc
md"Plot configuration for Ic, Vs and E₀"

# ╔═╡ aa4bfc30-6e61-48b9-805e-7d7088946de7
begin
	
	#Figure
	figIc = Figure(size = (600,600))
	
	#Axes
	Axis(figIc[1,1], xticks = (0:0.5:3.0), yticks = (ytick_min:0), ylabel = "Elevation (m)", xlabel = "Ic", yticklabelsize = axislabelsize, xticklabelsize = axislabelsize, limits = mylimits)
	Axis(figIc[1,2], xticks = (0:50:300), yticks = (ytick_min:0),  xlabel = "Vs (m/s)", yticklabelsize = axislabelsize, xticklabelsize = axislabelsize, limits = mylimits)
	Axis(figIc[1,3], xticks = (0:100:maximum(E0)), yticks = (ytick_min:0), xlabel = "E₀ (MPa)", yticklabelsize = axislabelsize, xticklabelsize = axislabelsize, limits = mylimits)
	
	#Plots
	lines!(figIc[1,1], Ic, -depth_m, color=RGBf(0.1, 0.1, 0.1))
	scatter!(figIc[1,1], get_soil_type_CPT2012(Ic), -depth_m, markersize = markersize, label = "Soil type for CPT2012")
	lines!(figIc[1,2], Vs, -depth_m, color=RGBf(0.1, 0.1, 0.1))
	lines!(figIc[1,3], E0, -depth_m, color=RGBf(0.1, 0.1, 0.1))
	lines!(figIc[1,3], fun_E0_linear(depth_m), -depth_m, linestyle = :dash)
	

		
end;

# ╔═╡ 8fa56b53-7b94-40a2-851a-63eb4014348c
	figIc

# ╔═╡ Cell order:
# ╠═64afe942-879b-4aea-84ae-cf1fcb766689
# ╠═ca156a62-ce63-11ef-32b1-bd2ade8dcf4c
# ╠═d1b05e7f-a074-43d9-be90-cd8381ea0c39
# ╟─7427ecd6-0dda-4284-948e-e5902c63daa1
# ╟─38c42433-425f-4f38-a99a-5ed2e4977ac9
# ╟─7f303256-ef8d-46e1-b21b-f02d9a4ec0fb
# ╟─387f30fb-3b91-405d-bb11-436f3e45f188
# ╟─ab327aef-981f-4898-be56-bfac4852ebd8
# ╟─7b404570-cebd-44af-8b87-7179cc2fb1f3
# ╟─d632faf0-75ee-489b-81b4-9b8575a083e2
# ╠═1e50db38-e0ad-4f88-9baa-8d708d577e5f
# ╟─ebe80e52-af49-4475-a32f-d61b36ec8716
# ╠═bf2c5ffc-d7d1-424d-9a4e-0f78b0ae5ebd
# ╟─8a3d5f67-4e4f-4c33-904c-5bde4a9cc3ef
# ╟─2a956b1c-be86-4218-a8d4-e6baafe048dd
# ╟─1b9997e7-a282-453d-b5a7-a4705bc4ad84
# ╠═8f3d0cc8-6466-43dd-82c8-4d09348e7fc5
# ╟─ac1b64ea-b0cb-4d81-8aa9-214c3bd3d42f
# ╟─b8301989-928b-409f-95bb-ca15fb5a2431
# ╠═db43b118-16b9-4d66-be86-a129fc286107
# ╟─8435f38b-91a6-47e1-b596-2296e16070f0
# ╠═88c706e2-7e38-4f3f-a257-16d522f964e7
# ╟─13a9ad72-1687-4b6c-998e-035be8e93918
# ╠═eb7b6073-d5b3-4bde-b390-5cc76b7d2cff
# ╟─7e9b0f08-1304-47f6-b981-9a4c001df324
# ╠═2df2143c-0e08-4bb6-b056-d56cdd4eb21a
# ╟─6edf2a4b-4672-437e-b8da-d0d91c9f924f
# ╟─a184c4eb-cb0e-4a99-8568-9b796ad0ab57
# ╠═2c47f0aa-22ea-4723-8c85-cadba3443472
# ╟─805713cf-e4fa-44cd-b694-5ebe05d7637f
# ╠═acb1bb78-6ff9-4b9a-85bd-b56c726635a1
# ╟─d8728f41-547d-4e4e-a1d0-689a53ac7d28
# ╟─46790341-9cae-455b-b0b1-231c10f5c6e1
# ╠═a8117f4a-a2fb-4a5f-bce2-0f88819a22a8
# ╟─7e38b7c1-5325-4c94-96c8-b06a2fc2c561
# ╠═715576a6-6f90-43c6-9e81-6af467429ae6
# ╟─ec1769aa-2424-4aed-8696-e157fab03673
# ╠═779a21a3-82b0-4d22-95dd-d1f896b3a272
# ╟─a3337224-351b-4776-9226-36ece5318f56
# ╟─ece64ce0-f597-4846-87d5-4146b7277260
# ╟─811a0b94-aaf9-4467-89e1-2723651dfa1c
# ╟─29115972-de54-48b7-b91e-3b747f249f4a
# ╠═679e440d-2e39-4cb7-a1c9-9c9ceeac9f6d
# ╟─e159c66d-36c2-409a-b287-f07e381f8c6a
# ╠═171bdb61-6e79-408a-8835-1706a415e688
# ╟─c166a7d0-39b8-4f03-bab7-984722cdfda2
# ╟─d6178e74-b9cb-4587-bc58-6bc1f42f9e2e
# ╟─12a5793a-30de-4f49-a95a-81bb0e4f48a3
# ╠═3d0b9618-f5ba-4b2e-a559-f5862d05b485
# ╟─d1024dbc-436e-44ed-96cc-1e030547575c
# ╠═5f5e2b35-629c-4a11-89a9-4b129894053b
# ╟─2c996085-e5b6-4d63-99db-6a47af2cd13b
# ╟─3968f88c-d0ea-4b6e-aea2-407584ad917e
# ╟─f7bb1887-3a38-4752-9ee2-769734bed772
# ╠═8414a932-8391-4e47-9912-f7039f1ccdb3
# ╟─7e99bf7f-4b8d-41a2-9319-887244f75a9c
# ╠═8f23be3d-ba6a-4340-8169-7881eaf7ddea
# ╟─93606b0b-7294-424a-8edc-3e7901fae6fd
# ╠═075ec3aa-666a-43d5-a551-10e6b9a10d3a
# ╟─f5546bf5-78e7-413e-b9f3-9bb6005438b6
# ╠═81a00b17-402d-4a20-a3ab-38b6eef32f6f
# ╟─f32428f3-b5c7-430f-bb2b-f9be098fcf53
# ╟─c10cfbb3-1f57-4884-a2bc-f2e5827e0b2e
# ╠═201ebe2b-fdaa-430d-b1d4-090f4c39a9b7
# ╟─3de17617-19d6-4550-8b4b-7bd2d38e3aee
# ╟─a1795fc2-61b1-4ccc-a39a-f3a18a5b38dd
# ╠═8fa56b53-7b94-40a2-851a-63eb4014348c
# ╟─97ed1519-b12e-4d72-980b-b5fca12d5cdb
# ╟─b4291248-c6f6-49ae-b702-b64169a96b3a
# ╟─5dca4557-ef1f-4737-a66a-c080decb6460
# ╟─e5a03703-8dfb-413a-8df9-67896cf9b2f4
# ╠═00a15171-61c0-4ed8-8db6-203cc288fba6
# ╟─3cb2c6ba-a463-4e58-9d3b-3d8577363784
# ╟─433be51b-8738-42d2-92c0-0f62344aab3f
# ╟─d980942c-e160-472b-98eb-a5a939e04d1a
# ╠═d8c6bcc9-355d-46ca-b13a-e130e8895ede
# ╠═fd8d59af-d979-4bcf-8aef-7ed338f95e71
# ╠═392f8de0-ae6b-40fe-9b77-cd4603831eec
# ╠═48be6b46-9061-483e-a574-4870d9f33669
# ╟─8cd9047c-3796-46f0-a27f-d98ef86888a6
# ╟─a30fb07b-ec4a-4e41-be61-70718894516a
# ╠═234d9744-0550-4485-82d1-b6a74ba6f1be
# ╠═673a87a4-569d-4c33-ae3e-98debc89df73
# ╠═d50e9d4c-2962-4735-bea8-4bd20d120f32
# ╟─35138641-da56-41d2-a266-3e0de45a832d
# ╠═f40fc671-1cd6-4b31-a1b3-a5da433f8805
# ╟─03e1dd77-669d-46b2-bff0-8e23fef7cf39
# ╠═380622fb-3e25-4259-9ae0-45af9fd88d0c
# ╠═dc029c10-1afb-4401-85c9-732aa7962510
# ╟─3a4fc7e3-cf0c-4ada-ae84-17de8ac74f52
# ╟─b518a716-062b-4e91-bff8-cbd51e09be1d
# ╟─15715b84-c935-4903-ae55-dc115c530dab
# ╟─7d5bfca5-5ebd-4e5a-8d7f-0e9c1684a922
# ╟─64013deb-af90-465d-a8f5-581e0f6baab6
# ╠═0b104fff-80fe-43e9-be05-99521aeb0535
# ╟─fe095cea-69a2-46a1-823c-dc1bb7613469
# ╟─efe4f0cc-b42b-4898-9103-0ca91747d1c8
# ╠═e5af32ab-6eb3-46d5-8a13-6d6efb8c8e9d
# ╠═4f9e4011-e02e-4fb7-a76f-d4ac0ca3ad43
# ╠═c232752c-4f7e-4b32-92f0-e793f4838975
# ╠═741eb051-a7cc-456a-b992-044d0acd460a
# ╠═8f9a73f7-d3ab-4f05-b5cc-7daeab7afd66
# ╠═4e9ccc3b-e977-4f15-b0c5-7c5daca84901
# ╠═c46f48bb-e03b-4b90-95cb-513039d5416f
# ╟─fef5eba2-7457-4bdf-a889-3650ba38335f
# ╟─c58cb696-e3f6-441d-81f1-5e4915db081d
# ╠═723a758c-941e-450f-ac2c-60b6a7e140db
# ╟─9ea12381-d2e4-442a-bb67-95fdf5f67437
# ╠═97e12fce-2e96-49fd-8476-0ea178de8c97
# ╟─8335941c-dc8a-4dd9-a558-5314c722c523
# ╠═338baa40-0742-4e48-9aca-6e7489c45e1b
# ╟─53f5866b-eecc-423a-9f64-b22a0bfc2de9
# ╠═bb86eff9-f3ea-4ade-ae3c-5feeac9d8e34
# ╟─6c5bdbe9-3253-4f6e-87aa-d62ba21016d2
# ╟─f2a332dd-c5fa-4a89-a2d7-89818f3f23b3
# ╟─b1b1c43f-cc35-4c18-a4ff-5ea0fb8996e6
# ╠═789627e2-0e67-4692-b750-a155b9f16790
# ╟─4e9ef2a1-e2a0-4209-8eb5-6212fb4fdb5f
# ╟─c0be21b4-a70a-4dc5-b2b9-dba8ff243a99
# ╟─dc5cb03e-2fda-4f0c-a580-0deada74707d
# ╠═650d88ea-5bdd-4fb1-b734-9184774681aa
# ╟─2fefead9-8aab-4044-92bb-2892b15836ca
# ╠═ff57b7dd-a174-4f15-9cd9-1bb08087c083
# ╟─53ef7ea8-93f5-466d-a7ae-66819ac8c1dc
# ╠═aa4bfc30-6e61-48b9-805e-7d7088946de7
