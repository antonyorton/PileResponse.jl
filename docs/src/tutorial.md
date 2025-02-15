```@meta
EditURL = "tutorial.jl"
```

# [Tutorial](@id tutorial)

Note: First add the `PileResponse` package if you have not already.

`Pkg> add "https://github.com/antonyorton/PileResponse.jl"`

````@example tutorial
import PileResponse as prs
using CairoMakie
CairoMakie.activate!(type="svg");
nothing #hide
````

## Read CPT data from file

````@example tutorial
data = prs.read_delimited_text_file("../../data/example_cpt_data.csv", delim=',');
nothing #hide
````

Extract the column names and then assign to variables

````@example tutorial
depth_col, qc_col, fs_col, u2_col = prs.find_cpt_column_names(collect(keys(data)));
depth_m = data[depth_col];
qc_MPa = data[qc_col];
fs_MPa = data[fs_col] * 0.001;
u2_MPa = data[u2_col] * 0.001;
nothing #hide
````

### Plot the CPT data

````@example tutorial
f = Figure(size=(900, 750));
ax1 = Axis(f[1, 1], title="qc (MPa)", xlabel="qc (MPa)", ylabel="Elevation (m)");
ax2 = Axis(f[1, 2], title="fs (kPa)", xlabel="fs (kPa)");
ax3 = Axis(f[1, 3], title="u2 (kPa)", xlabel="u2 (kPa)");
lines!(f[1, 1], qc_MPa, -depth_m);
lines!(f[1, 2], fs_MPa, -depth_m);
lines!(f[1, 3], qc_MPa, -depth_m);
f
````

## Soil properties and depth to groundwater

Soil properties

````@example tutorial
gamma_soil = 20.5; #kN/m3
Poisson_ratio = 0.3;
nothing #hide
````

Groundwater depth

````@example tutorial
gw_depth = 3.0; # m below ground
nothing #hide
````

## Derive ``I_{c}``, ``V_{s}`` and ``E_{0}``

````@example tutorial
Ic = prs.get_Ic(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73);
Vs = prs.get_Vs(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73);
E0 = prs.get_E0(Vs, gamma=gamma_soil, ν=Poisson_ratio);
nothing #hide
````

### Plot the derived values from the CPT data

````@example tutorial
f = Figure(size=(900, 750));
ax1 = Axis(f[1, 1], title="Ic", xlabel="Ic", ylabel="Elevation (m)");
ax2 = Axis(f[1, 2], title="Vs (m/s)", xlabel="Vs (m/s)");
ax3 = Axis(f[1, 3], title="E₀ (MPa)", xlabel="E₀ (MPa)");
lines!(f[1, 1], Ic, -depth_m, color=RGBf(0.1, 0.5, 0.1));
lines!(f[1, 2], Vs, -depth_m, color=RGBf(0.1, 0.5, 0.1));
lines!(f[1, 3], E0, -depth_m, color=RGBf(0.1, 0.5, 0.1));
f
````

Some text below the figure.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

