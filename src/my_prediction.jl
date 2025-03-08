
import PileResponse as prs
using PrettyTables
using CairoMakie
using Printf
CairoMakie.activate!(type="svg");

"""This file provides an example prediction"""

cpt_filename = "./data/example_cpt_data.csv"
test_id = cpt_filename[findlast("/", cpt_filename)[1]+1:end][1:end-4] #extract the file name minus the extension

# load_test_results_file = "./data/load_test_results.csv"; #include load test results if available

gamma_soil = 20.5; #Soil unit weight (kN/m3)
Poisson_ratio = 0.3;
gw_depth = 3.00; # m below ground
pile_type = "Driven pile - steel closed ended"
pile_length = 17.0; # metres
pile_diameter = 0.457; # metres (avg of equal area and equal perimeter circle for the square pile)
Epile_MPa = 27400 #MPa

data = prs.read_delimited_text_file(cpt_filename, delim=',');

depth_col, qc_col, fs_col, u2_col = prs.find_cpt_column_names(collect(keys(data)));
depth_m = data[depth_col];
qc_MPa = data[qc_col];
fs_MPa = data[fs_col] * 0.001;
u2_MPa = data[u2_col] * 0.001;


# Get Ic, Vs, E0
Ic = prs.get_Ic(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73);
Vs = prs.get_Vs(depth_m, qc_MPa, fs_MPa, u2_MPa, gw_depth, gamma=gamma_soil, a=0.73);
E0 = prs.get_E0(Vs, gamma=gamma_soil, ν=Poisson_ratio);

# Least squares approx to E0
fun_E0_linear = prs.get_least_squares_interpolator(depth_m, E0);

@printf "The pile type is: %s\n" pile_type

########
## Pile ultimate load

#Get soil types
soil_type_CPT2012 = prs.get_soil_type_CPT2012(Ic);

#Ultimate shaft load
fshaft_MPa = prs.get_ultimate_shaft_resistance(qc_MPa, Ic, pile_type, factor=1.0);
ult_shaft_MN = prs.get_ultimate_shaft_load(depth_m, fshaft_MPa, pile_diameter, pile_length);

@printf "The ultimate shaft load is %5.3f MN\n" ult_shaft_MN

#Ultimate base load
qc_avg_base = prs.get_average_qc_at_pile_base(depth_m, qc_MPa, pile_length, pile_diameter, clip_to_30pct=true);

kc_at_base = prs.get_kc_base_CPT2012()[pile_type][soil_type_CPT2012[argmin(abs.(depth_m .- pile_length))]][1];

fb_MPa = kc_at_base * qc_avg_base;
ult_base_MN = pi * pile_diameter^2 / 4 * fb_MPa;

@printf "qc_avg_base = %5.3f MPa\n" qc_avg_base
@printf "kc_at_base = %5.2f \n" kc_at_base
@printf "The ultimate base load is %5.3f MN\n" ult_base_MN

#Pile ultimate load
pile_ult_load = ult_base_MN + ult_shaft_MN;

@printf "The pile ultimate load is %5.3f MN\n" pile_ult_load


########
### Load displacement response

# Stiffness parameters
E_L = fun_E0_linear(pile_length)
E_Lon2 = fun_E0_linear(pile_length / 2)
k0 = prs.get_initial_pile_head_stiffness(pile_length, pile_diameter, Epile_MPa, E_L, E_Lon2, ν=Poisson_ratio)

@printf "E₀ at the base of the shaft = %4.0f MPa\n" E_L
@printf "E₀ at the midpoint of the shaft = %4.0f MPa\n" E_Lon2
@printf "k₀ = %4.0f MN/m\n" k0

# Loads and displacement
pile_head_loads = 0.01:0.001:0.9*pile_ult_load;
displacement = prs.get_pile_head_displacement(k0, pile_head_loads, pile_ult_load);

# # Pile capacity
# index_capacity = argmin(abs.(displacement .- allowable_pile_head_settlement_m))
# allowable_pile_head_settlement_m = 0.03;
# pile_capacity_MN = pile_head_loads[index_capacity];
# disp_at_capacity = displacement[index_capacity]

# @printf "The pile capacity is %5.3f MN\n" pile_capacity_MN
# @printf "The displacement at capacity is %5.3f m\n" disp_at_capacity


# Load test results (for comparison)
# load_test = prs.read_delimited_text_file("./data/load_test_results.csv");


#Plot CPT raw data
fig1 = Figure(size=(800, 800));
ax1 = Axis(fig1[1, 1], title=string(test_id, " - qc (MPa)"), xlabel="qc (MPa)", ylabel="Elevation (m)", yticks=(-round(depth_m[end] + 1):0));
ax2 = Axis(fig1[1, 2], title="fs (kPa)", xlabel="fs (kPa)", yticks=(-round(depth_m[end] + 1):0));
ax3 = Axis(fig1[1, 3], title="u2 (kPa)", xlabel="u2 (kPa)", yticks=(-round(depth_m[end] + 1):0));
lines!(fig1[1, 1], qc_MPa, -depth_m);
lines!(fig1[1, 2], 1000 * fs_MPa, -depth_m);
lines!(fig1[1, 3], 1000 * u2_MPa, -depth_m);
save("figure1_raw_data.pdf", fig1, pdf_version="1.4")


#Plot CPT derived data
fig2 = Figure(size=(800, 800));
ax1 = Axis(fig2[1, 1], title=string(test_id, " - Ic"), xlabel="Ic", ylabel="Elevation (m)", yticks=(-round(depth_m[end] + 1):0));
ax2 = Axis(fig2[1, 2], title="Vs (m/s)", xlabel="Vs (m/s)", yticks=(-round(depth_m[end] + 1):0));
ax3 = Axis(fig2[1, 3], title="E₀ (MPa)", xlabel="E₀ (MPa)", yticks=(-round(depth_m[end] + 1):0));
lines!(fig2[1, 1], Ic, -depth_m);
lines!(fig2[1, 2], Vs, -depth_m);
lines!(fig2[1, 3], E0, -depth_m);
lines!(fig2[1, 3], fun_E0_linear.(depth_m), -depth_m, linestyle=:dash, color=:black);
save("figure2_derived_data.pdf", fig2, pdf_version="1.4")


#Plot the load displacement response to capacity
fig3 = Figure(size=(800, 500));
ax1 = Axis(fig3[1, 1], title=string(test_id, " - Pile head response"), xlabel="Displacement (m)", ylabel="Load (MN)");
lines!(fig3[1, 1], displacement, pile_head_loads, label="My prediction")
# lines!(fig3[1, 1], load_test["displacement_m"], load_test["load_MN"], color=:red, label="Load test record")
axislegend("Legend", position=:lt)
save("figure3_load_response.pdf", fig3, pdf_version="1.4")