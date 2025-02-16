```@meta
EditURL = "tutorial.jl"
```

# [Tutorial](@id tutorial)

Note: Remember to add the `PileResponse` package if you have not already.
 - `Pkg> add "https://github.com/antonyorton/PileResponse.jl"`

````@example tutorial
push!(LOAD_PATH, "../../src/") #This line can be ignored if the package is added as above
import PileResponse as prs
using CairoMakie
CairoMakie.activate!(type="svg");
nothing #hide
````

## Read CPT data
The first thing we need to do is to read CPT data from a text file such as a .csv file. It is a requirement that the file has columns for:
- depth (m)
- qc (MPa)
- fs (kPa)
- u2 (kPa)

The column names should be sensible or the program will not be able to guess them.
You will also need to ensure the correct relative path to your file and the correct delimeter, for example comma (',') or tab ('\t').

````@example tutorial
data = prs.read_delimited_text_file("../../data/example_cpt_data.csv", delim=',');
nothing #hide
````

We now extract the column names and assign the data to variables.

````@example tutorial
depth_col, qc_col, fs_col, u2_col = prs.find_cpt_column_names(collect(keys(data)));
depth_m = data[depth_col];
qc_MPa = data[qc_col];
fs_MPa = data[fs_col] * 0.001;
u2_MPa = data[u2_col] * 0.001;
nothing #hide
````

### Plot the CPT data
It is recommended to plot the data and check that it is as expected by comparing to pdf reports if available.

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
We assign some properties here which are assumed to be constant over the soil profile. The loss in accuracy from assuming a constant soil unit weight is likely to be acceptable for most situations.

Soil properties

````@example tutorial
gamma_soil = 20.5; #Soil unit weight (kN/m3)
Poisson_ratio = 0.3;
nothing #hide
````

Groundwater depth

````@example tutorial
gw_depth = 3.0; # m below ground
nothing #hide
````

## Derive ``I_{c}``, ``V_{s}`` and ``E_{0}``
From the correlations published by Roberston and Cabal (2022), we derive\
\
The soil behaviour type
  - ``I_{c} = \sqrt{(3.47 - log_{10}(Q_{t}))^{2} + (log_{10}(F_{r}) + 1.22)^{2}}`` , and
The shear wave velocity
  - ``\sqrt{(\alpha_{vs}\cdot \large\frac{q_{n}}{\small{0.101}})}``, where ``\alpha_{vs} = 10^{(0.55 Ic + 1.68)}``
We then use the following two relationships to derive the small strain elastic modulus $E_{0}$:
  - ``G_{0} = \frac{\gamma}{9.81} \cdot V_{s}^2``

  - ``E_{0} = 2(1 + ν)G_{0}``

````@example tutorial
Ic = prs.get_Ic(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73);
Vs = prs.get_Vs(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73);
E0 = prs.get_E0(Vs, gamma=gamma_soil, ν=Poisson_ratio);
nothing #hide
````

### Linear approximation to ``E_{0}`` with depth
For the assessment of the intial pile head stiffness, we need to model the soil profile as having a linearly increasing elastic modulus with depth. This is common in many situations.

We fit a line to the $E_{0}$ data using a least squares fit.

````@example tutorial
fun_E0_linear = prs.get_least_squares_interpolator(depth_m, E0);
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
lines!(f[1, 3], fun_E0_linear(depth_m), -depth_m, color=RGBf(0.1, 0.5, 0.1), linestyle=:dash);
f
````

## Pile load test details
For this tutorial, we consider a pile load test which was carried out at Mobile, Alabama. The test comprised a 457 mm diameter, 9.5 mm wall, closed-toe, steel-pipe pile driven to 17.0 m below ground and grouted after driving. The pile stick up was approximately 0.9 m above ground.\
\
Further details of the test can be found at the website of the test organiser [www.fellenius.net](https://www.fellenius.net/)

### Pile properties
We assign a pile type from the list of available pile types.

````@example tutorial
prs.list_available_pile_types()
````

````@example tutorial
pile_type = "Driven pile - steel closed ended";
nothing #hide
````

````@example tutorial
pile_toe_depth = 17.0; # metres
nothing #hide
````

````@example tutorial
pile_diameter = 0.457; # metres
nothing #hide
````

The pile elastic modulus `Epile_MPa` for this case is caculated assuming a steel shell, grout filled circular pile

````@example tutorial
odiam = 0.457;
idiam = 0.457 - 0.0095;
Aouter = pi * odiam^2 / 4;
Ainner = pi * idiam^2 / 4;
Eouter = 200 # GPa;
Einner = 20  # GPa;
Epile_GPa = ((Aouter - Ainner) * Eouter + Ainner * Einner) / Aouter;
nothing #hide
````

````@example tutorial
Epile_MPa = round(Int64, 1000 * Epile_GPa);
print("The elastic modulus of the pile is ", Epile_MPa, " (MPa)")
````

## Pile ultimate load
In this section we calculate the pile ultimate load following the approach of Frank (2017).\
\

**Soil types**\
The calculation of pile ultimate resistance relies on several factors which require a soil type.
Using the assessed ``I_{c}`` value, the soil types required for the assessment are defined in this notebook as follows:
- ``I_{c} > 2.60`` : Silt and clay (soil type 1).
- ``2.05 < I_{c} \leq 2.60`` : Intermediate soil (soil type 2).
- ``I_{c} \leq 2.05`` : Sand (soil type 3).

````@example tutorial
soil_type_CPT2012 = prs.get_soil_type_CPT2012(Ic);
nothing #hide
````

### Pile ultimate shaft load
The ultimate shaft resistance for each node along the pile shaft is obtained as follows:
- ``f_{s} = \alpha\cdot f_{sol}``, with ``f_{sol} \leq f_{smax}``\
\
where:
- ``f_{sol}`` is the unfactored shaft resistance dependent on ``q_{c}`` and soil type.
- ``\alpha`` is a factor dependent on soil type and pile type.\
- ``f_{smax}`` is a the limiting shaft resistance dependent on soil type and pile type.

First, calculate the ultimate resistance (MPa) for each node along the pile

````@example tutorial
fshaft_MPa = prs.get_ultimate_shaft_resistance(qc_MPa, Ic, pile_type, factor=1.0);
nothing #hide
````

Then we can calculate the ultimate shaft load (MN) for the whole pile

````@example tutorial
ult_shaft_MN = prs.get_ultimate_shaft_load(depth_m, fshaft_MPa, pile_diameter, pile_toe_depth);
print("The ultimate shaft load is ", round(ult_shaft_MN, digits=3), " (MN)")
````

### Pile ultimate base load
The ultimate base load is obtained as follows:\
\
First, we obtain the ultimate base resistance (MPa):
- ``f_{b} = k_{c}\cdot q_{ca}``\
where:
- ``k_{c}`` is a factor dependent on soil type and pile class, and;\
- ``q_{ca}`` is the equavalent average cone resistance at the pile base.

````@example tutorial
qc_avg_base = prs.get_average_qc_at_pile_base(depth_m, qc_MPa, pile_toe_depth, pile_diameter, clip_to_30pct=false);
nothing #hide
````

````@example tutorial
kc_at_base = prs.get_kc_base_CPT2012()[pile_type][soil_type_CPT2012[depth_m.==pile_toe_depth]][1];
nothing #hide
````

We then caculate the ultimate base load.

````@example tutorial
fb_MPa = kc_at_base * qc_avg_base;
ult_base_MN = pi * pile_diameter^2 / 4 * fb_MPa;
nothing #hide
````

Display the results

````@example tutorial
print("The average qc within 1.5 diameters from the pile base is ", round(qc_avg_base, digits=2), " (MPa)", "\n")
print("The factor kc is ", kc_at_base, "\n")
print("The ultimate base load is ", round(ult_base_MN, digits=3), " (MN)", "\n")
````

### Predicted pile ultimate load

````@example tutorial
print("The predicted pile ultimate load is ", round(ult_base_MN + ult_shaft_MN, digits=3), "(MPa)")
````

\
Further sections to come soon ....

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

