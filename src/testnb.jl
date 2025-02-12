### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ d27efcff-916b-44f8-ad8f-4d8d08ae780e
# ╠═╡ show_logs = false
begin
    import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.add(url = "https://github.com/antonyorton/PileResponse.jl")
    import PileResponse as prs
end;

# ╔═╡ 07a91a4c-f877-4d13-a9b6-61938a4848eb
begin
	using CairoMakie
	CairoMakie.activate!(type = "svg");
end;

# ╔═╡ dfc0d050-e922-11ef-06ae-b1c7dbc6526f
md"""## Pile response from CPT
	**Author:** Antony Orton"""

# ╔═╡ 03d5a5fa-372f-4cd4-baed-d65594db62cb
md"""---------------------
Let's read a CPT data file
"""

# ╔═╡ 4d319fda-a82a-483b-a780-1874cf00c181
cptfilename = "../data/example_cpt_data.csv"

# ╔═╡ 02af09f4-db32-408f-9283-6fc77d12d792
data = prs.read_delimited_text_file(cptfilename, delim=',');

# ╔═╡ 24cbc2c5-ce1e-442a-9813-2e78767d3ed1
md"""-------------------------
Now we extract the column names and assign data to variables
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
Plot the CPT data"""

# ╔═╡ d73995c2-1e61-4f19-a02c-3ed89ad13fbb
md"""------------
Get Ic, Vs and E₀
"""

# ╔═╡ b7650638-fd37-40e8-aa82-59bdfb6da2a1
md"""
Input some soil properties
"""

# ╔═╡ 69c58569-4fe9-403a-a67e-817fc2911cc7
begin
	gamma_soil = 20.5; # kN/m3
	Poisson_ratio = 0.3;
	groundwater_depth = 3.0; # m below ground
end

# ╔═╡ 3c438add-3c28-41b1-824b-4e15c8b70f7d
begin
	Ic = prs.get_Ic(depth_m, qc_MPa, u2_MPa, groundwater_depth, gamma = gamma_soil, a = 0.73)
end

# ╔═╡ 013214c8-4475-4b92-8149-31f6dee0b7dd


# ╔═╡ 74da7044-609d-48a9-8ec3-63624d5ad87d


# ╔═╡ ca480ba8-b1e7-434c-aaa7-d33e3bb230ad


# ╔═╡ 4c8139ac-c625-4168-b86d-e1c3a3a0561e


# ╔═╡ 40c4b45a-d2ed-422b-b76e-2f3aa83b624b


# ╔═╡ c6730ab0-4669-4ff4-9658-df539fe4169d
begin
	#Common plot parameters
	ytick_min = -floor(round(depth_m[end], digits = 2)) - 1.0; #ymin
	mylimits = (0.0, nothing, ytick_min, 0.001) #(xmin,xmax,ymin,ymax)
	markersize = 4.0;
	rasterize = 3;
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

# ╔═╡ Cell order:
# ╟─dfc0d050-e922-11ef-06ae-b1c7dbc6526f
# ╠═d27efcff-916b-44f8-ad8f-4d8d08ae780e
# ╠═07a91a4c-f877-4d13-a9b6-61938a4848eb
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
# ╠═b7650638-fd37-40e8-aa82-59bdfb6da2a1
# ╠═69c58569-4fe9-403a-a67e-817fc2911cc7
# ╠═3c438add-3c28-41b1-824b-4e15c8b70f7d
# ╠═013214c8-4475-4b92-8149-31f6dee0b7dd
# ╠═74da7044-609d-48a9-8ec3-63624d5ad87d
# ╠═ca480ba8-b1e7-434c-aaa7-d33e3bb230ad
# ╠═4c8139ac-c625-4168-b86d-e1c3a3a0561e
# ╠═40c4b45a-d2ed-422b-b76e-2f3aa83b624b
# ╠═c6730ab0-4669-4ff4-9658-df539fe4169d
# ╠═48366eb0-8138-4719-a0dd-e9370d9cb806
