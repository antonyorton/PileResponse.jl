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

# ╔═╡ d27efcff-916b-44f8-ad8f-4d8d08ae780e
# ╠═╡ show_logs = false
begin
    import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.add(url = "https://github.com/antonyorton/PileResponse.jl")
    Pkg.add("PlutoUI")
	Pkg.add("CairoMakie")
	import PileResponse as prs
	using PlutoUI
	using CairoMakie
	CairoMakie.activate!(type = "svg");
end;

# ╔═╡ dfc0d050-e922-11ef-06ae-b1c7dbc6526f
md"""# Pile response from CPT
	**Author:** Antony Orton"""

# ╔═╡ 03d5a5fa-372f-4cd4-baed-d65594db62cb
md"""---------------------
### Read a CPT data file
"""

# ╔═╡ 4d319fda-a82a-483b-a780-1874cf00c181
cptfilename = "../data/example_cpt_data.csv"

# ╔═╡ 02af09f4-db32-408f-9283-6fc77d12d792
data = prs.read_delimited_text_file(cptfilename, delim=',');

# ╔═╡ 24cbc2c5-ce1e-442a-9813-2e78767d3ed1
md"""-------------------------
### Extract the column names and assign data to variables
"""

# ╔═╡ 5e0cacc3-b7ff-4398-bd1b-89e071127053
columns = collect(keys(data));

# ╔═╡ 1ba92e87-98a1-4207-b97f-0835d1c02e2d
begin
	depth_col = [item for item in columns if occursin("depth",lowercase(item))][1]
	qc_col = [item for item in columns if occursin("qc",lowercase(item))][1]
	fs_col = [item for item in columns if occursin("fs",lowercase(item))][1]
	u2_col = [item for item in columns if occursin("u2",lowercase(item))][1]
end;

# ╔═╡ eb61fe89-d2e9-4445-a529-86448e9b6dc5
[depth_col, qc_col, fs_col, u2_col]

# ╔═╡ 8470d353-7978-4769-ba19-3e6cc2648fac
md"""Assign to variables"""

# ╔═╡ c629a73a-ae6b-4702-a5c8-325add91b85e
begin
	depth_m = data[depth_col]
	qc_MPa = data[qc_col]
	fs_MPa = data[fs_col] * 0.001
	u2_MPa = data[u2_col] * 0.001
end;

# ╔═╡ 1012192d-2fd6-4747-bc61-f11d45d84b17
md"""---------------
### Plot the CPT data"""

# ╔═╡ d73995c2-1e61-4f19-a02c-3ed89ad13fbb
md"""------------
### Derive Ic, Vs and E₀
"""

# ╔═╡ b7650638-fd37-40e8-aa82-59bdfb6da2a1
md"""
Input soil properties and depth to groundwater
"""

# ╔═╡ 69c58569-4fe9-403a-a67e-817fc2911cc7
begin
	gamma_soil = 20.5; # kN/m3
	Poisson_ratio = 0.3;
	groundwater_depth = 3.0; # m below ground
end;

# ╔═╡ e6d04d0a-eb5e-481c-a8a6-b35ede6aca22
md"""Calculate Ic, Vs and E₀
"""

# ╔═╡ 3c438add-3c28-41b1-824b-4e15c8b70f7d
begin
	Ic = prs.get_Ic(depth_m, qc_MPa, fs_MPa, u2_MPa, groundwater_depth, gamma = gamma_soil, a = 0.73)
	Vs = prs.get_Vs(depth_m, qc_MPa, fs_MPa, u2_MPa, groundwater_depth, gamma = gamma_soil, a = 0.73)
	E0 = prs.get_E0(Vs, gamma = gamma_soil, ν = 0.3)
end;

# ╔═╡ 013214c8-4475-4b92-8149-31f6dee0b7dd
md"""Get a least squares fit to E₀"""

# ╔═╡ 74da7044-609d-48a9-8ec3-63624d5ad87d
fun_E0_linear = prs.get_least_squares_interpolator(depth_m, E0);

# ╔═╡ ca480ba8-b1e7-434c-aaa7-d33e3bb230ad
md"""Plot the derived results"""

# ╔═╡ 40c4b45a-d2ed-422b-b76e-2f3aa83b624b
md"""-----------------
### Pile ultimate load
"""

# ╔═╡ 1669ec82-23f0-48b2-85b1-3cadeb2c9f94
md"""Input the pile details"""

# ╔═╡ feb3b2e6-5804-4234-b59e-35ae9372a5bd
begin
	pile_toe_depth = 17.0 	# metres
	pile_diameter = 0.457 	# metres
	Epile_MPa = 27400 		# MPa
	allowable_pile_head_settlement_m = 0.060 	# metres
end;

# ╔═╡ 6122e0d1-78e7-432e-822f-45f67dda3583
md"Select the pile type"

# ╔═╡ 60a91466-170a-4d11-a412-26be9ba590d0
@bind pile_type Select(sort(collect(keys(prs.get_alpha_shaft_CPT2012()))))

# ╔═╡ a68d325a-bfe1-40dc-90a9-40eedb5f3701
md"""Get the soil types for the CPT2012 method (Frank, 2017)"""

# ╔═╡ 4ee315b6-edbe-4197-873d-c69f804870c2
soil_type_CPT2012 = prs.get_soil_type_CPT2012(Ic);

# ╔═╡ f3e79011-4d3d-4fe4-93b0-2b0bab5eea6f
md" #### Ultimate shaft load"

# ╔═╡ e94bd027-e536-44ff-9bc8-a22d7d40f4e0
md"Get the ultimate shaft resistance (MPa) for each node"

# ╔═╡ b6ae3d81-38cb-4a80-9219-3833d4463666
fshaft_MPa = prs.get_ultimate_shaft_resistance(qc_MPa, Ic, pile_type, factor = 1.0)

# ╔═╡ b95a595e-1c89-4a31-9738-1e4bc5ad4c76
md"Calulate the ultimate shaft load (MN) for the pile"

# ╔═╡ 5a322c8e-01e7-426f-93e3-ee98c6efb008
begin
	sum = 0.0; depth = 0.0; count = 1; #initial values
	while depth < pile_toe_depth
		#increment the sum for each shaft element
		sum += pi * pile_diameter * (depth_m[count + 1] - depth_m[count]) *
			0.5 * (fshaft_MPa[count + 1] + fshaft_MPa[count])
		#update the current depth and the count
		depth = depth_m[count + 1]
		count += 1
	end
	ult_shaft_MN = sum;
end;

# ╔═╡ a20911f8-4211-40ad-8600-4731bf4db194
md"The ultimate shaft load is **$(round(ult_shaft_MN, digits = 3)) (MN)**"

# ╔═╡ c0291985-31f2-42d7-b704-57a814bf79fd
md" #### Ultimate base load"

# ╔═╡ ed194582-d9dc-429e-8c90-51f5ff1dd70e


# ╔═╡ 9b13cd35-23c4-47fd-99c6-6bb0c1465dd4
md"Get the kc values for the selected pile type"

# ╔═╡ 17271f9f-2304-40f6-b18d-3a0e2f15a3fc
kc = prs.get_kc_base_CPT2012()[pile_type]

# ╔═╡ 15071d39-4f92-4150-9e19-89515174a9e8
md"Then for the soil type at the base of the pile"

# ╔═╡ f1b72e6b-35eb-4729-a337-10293db44411
kc_at_base = kc[soil_type_CPT2012[depth_m .== pile_toe_depth]][1]

# ╔═╡ 3a1328fb-4435-4589-b6e1-9117823bc92f
md"Get the average qc value within +/- 1.5 diameters of the pile base"

# ╔═╡ b40eb8ce-7524-46bd-988d-34251ef78d0c
qc_avg_base = prs.get_average_qc_at_pile_base(depth_m, qc_MPa, pile_toe_depth, pile_diameter, clip_to_30pct = false)

# ╔═╡ 57d3efc5-531b-43f1-a859-e6f11a0c266b
md"---------------------
### Plot configurations
"

# ╔═╡ c6730ab0-4669-4ff4-9658-df539fe4169d
begin
	#Common plot parameters
	ytick_min = -floor(round(depth_m[end], digits = 2)) - 1.0; #ymin
	mylimits = (0.0, nothing, ytick_min, 0.001) #(xmin,xmax,ymin,ymax)
	markersize = 3.0;
	rasterize = 1;
	axislabelsize = 10
	defaultlinecolour = color = RGBAf(0.2,0.6,0.2)
end;

# ╔═╡ 48366eb0-8138-4719-a0dd-e9370d9cb806
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
	# #smoothed data
	# lines!(fig_cpt[1,1], fun_qc(new_depth), -new_depth, rasterize = rasterize, color=RGBf(0.1, 0.1, 0.1));
	# lines!(fig_cpt[1,2], fun_fs(new_depth), -new_depth, rasterize = rasterize, color=RGBf(0.1, 0.1, 0.1));
	# lines!(fig_cpt[1,3], fun_u2(new_depth), -new_depth, rasterize = rasterize, color=RGBf(0.1, 0.1, 0.1));
	# save("figure.pdf", fig_cpt, pdf_version="1.4");
end;

# ╔═╡ b433cd7f-ace0-493e-bc2b-94ea7abf8804
fig_cpt

# ╔═╡ 5a24654d-fc65-429f-b6e7-eaf750f38356
begin
	
	#Figure
	figIc = Figure(size = (600,600))
	
	#Axes
	Axis(figIc[1,1], xticks = (0:0.5:3.0), yticks = (ytick_min:0), ylabel = "Elevation (m)", xlabel = "Ic", yticklabelsize = axislabelsize, xticklabelsize = axislabelsize, limits = mylimits)
	Axis(figIc[1,2], xticks = (0:50:300), yticks = (ytick_min:0),  xlabel = "Vs (m/s)", yticklabelsize = axislabelsize, xticklabelsize = axislabelsize, limits = mylimits)
	Axis(figIc[1,3], xticks = (0:100:maximum(E0)), yticks = (ytick_min:0), xlabel = "E₀ (MPa)", yticklabelsize = axislabelsize, xticklabelsize = axislabelsize, limits = mylimits)
	
	#Plots
	lines!(figIc[1,1], Ic, -depth_m, color=RGBf(0.1, 0.1, 0.1))
	scatter!(figIc[1,1], prs.get_soil_type_CPT2012(Ic), -depth_m, markersize = markersize, label = "Soil type for CPT2012")
	lines!(figIc[1,2], Vs, -depth_m, color=RGBf(0.1, 0.1, 0.1))
	lines!(figIc[1,3], E0, -depth_m, color=RGBf(0.1, 0.1, 0.1))
	lines!(figIc[1,3], fun_E0_linear(depth_m), -depth_m, linestyle = :dash)
	
end;

# ╔═╡ 4c8139ac-c625-4168-b86d-e1c3a3a0561e
figIc

# ╔═╡ Cell order:
# ╟─dfc0d050-e922-11ef-06ae-b1c7dbc6526f
# ╠═d27efcff-916b-44f8-ad8f-4d8d08ae780e
# ╟─03d5a5fa-372f-4cd4-baed-d65594db62cb
# ╠═4d319fda-a82a-483b-a780-1874cf00c181
# ╠═02af09f4-db32-408f-9283-6fc77d12d792
# ╟─24cbc2c5-ce1e-442a-9813-2e78767d3ed1
# ╠═5e0cacc3-b7ff-4398-bd1b-89e071127053
# ╠═1ba92e87-98a1-4207-b97f-0835d1c02e2d
# ╠═eb61fe89-d2e9-4445-a529-86448e9b6dc5
# ╟─8470d353-7978-4769-ba19-3e6cc2648fac
# ╠═c629a73a-ae6b-4702-a5c8-325add91b85e
# ╟─1012192d-2fd6-4747-bc61-f11d45d84b17
# ╠═b433cd7f-ace0-493e-bc2b-94ea7abf8804
# ╟─d73995c2-1e61-4f19-a02c-3ed89ad13fbb
# ╟─b7650638-fd37-40e8-aa82-59bdfb6da2a1
# ╠═69c58569-4fe9-403a-a67e-817fc2911cc7
# ╟─e6d04d0a-eb5e-481c-a8a6-b35ede6aca22
# ╠═3c438add-3c28-41b1-824b-4e15c8b70f7d
# ╟─013214c8-4475-4b92-8149-31f6dee0b7dd
# ╠═74da7044-609d-48a9-8ec3-63624d5ad87d
# ╟─ca480ba8-b1e7-434c-aaa7-d33e3bb230ad
# ╠═4c8139ac-c625-4168-b86d-e1c3a3a0561e
# ╟─40c4b45a-d2ed-422b-b76e-2f3aa83b624b
# ╟─1669ec82-23f0-48b2-85b1-3cadeb2c9f94
# ╠═feb3b2e6-5804-4234-b59e-35ae9372a5bd
# ╟─6122e0d1-78e7-432e-822f-45f67dda3583
# ╟─60a91466-170a-4d11-a412-26be9ba590d0
# ╟─a68d325a-bfe1-40dc-90a9-40eedb5f3701
# ╠═4ee315b6-edbe-4197-873d-c69f804870c2
# ╟─f3e79011-4d3d-4fe4-93b0-2b0bab5eea6f
# ╟─e94bd027-e536-44ff-9bc8-a22d7d40f4e0
# ╠═b6ae3d81-38cb-4a80-9219-3833d4463666
# ╟─b95a595e-1c89-4a31-9738-1e4bc5ad4c76
# ╠═5a322c8e-01e7-426f-93e3-ee98c6efb008
# ╟─a20911f8-4211-40ad-8600-4731bf4db194
# ╟─c0291985-31f2-42d7-b704-57a814bf79fd
# ╠═ed194582-d9dc-429e-8c90-51f5ff1dd70e
# ╟─9b13cd35-23c4-47fd-99c6-6bb0c1465dd4
# ╠═17271f9f-2304-40f6-b18d-3a0e2f15a3fc
# ╟─15071d39-4f92-4150-9e19-89515174a9e8
# ╟─f1b72e6b-35eb-4729-a337-10293db44411
# ╠═3a1328fb-4435-4589-b6e1-9117823bc92f
# ╠═b40eb8ce-7524-46bd-988d-34251ef78d0c
# ╟─57d3efc5-531b-43f1-a859-e6f11a0c266b
# ╠═c6730ab0-4669-4ff4-9658-df539fe4169d
# ╠═48366eb0-8138-4719-a0dd-e9370d9cb806
# ╠═5a24654d-fc65-429f-b6e7-eaf750f38356
