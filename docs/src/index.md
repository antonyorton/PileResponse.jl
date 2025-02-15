# [PileResponse.jl Documentation](@id PileResponse.jl)

```@contents
Pages = ["index.md", "tutorial.md", "api.md"]
Depth = 3
```

`Pkg> add "https://github.com/antonyorton/PileResponse.jl"`

PileResponse.jl is a module designed for the assessment of the load-deflection response of a pile head under static vertical load. The assessment depends on the pile type and the results of a single cone penetrometer test (CPT) conducted to at least the design toe level of the pile.\
\
The assessment of ultimate pile shaft and base capacities follows the approach outlined by Frank (2017). The assessment of the non linear load-deflection response is based on the initial pile head stiffness estimated by the Randolph and Wroth (1978) formula. The pile head stiffness is then factored down for increasing loads following the equation proposed by Mayne (2001).\
\
The interpretation of soil behaviour type ($I_{c}$) and shear wave velocity ($V_{s}$) from the CPT data follow correlations provided by Robertson and Cabal (2022).

**References**

Frank, R. (2017). Some aspects of research and practice for pile design in France. Innov. Infrastruct. Solutions. 2 (32) 1 - 15.\
\
Fahey, M. and Carter, J. P. (1993). A finite element study of the pressuremeter in sand using a nonlinear elastic plastic model. Can. Geot. Jnl. 30(2): 348 - 362.\
\
Randolph, M. F. and Wroth, C. P. (1978). Analysis of deformation of vertically loaded piles. Jnl. Geot. Eng. Divn., ASCE, 108 (GT12): 1465 - 1488.\
\
Mayne, P. W. (2001). Stress-strain-strength-flow parameters from enhanced in-situ tests. Proc. Int. Conf. In-Situ Measurement of Soil Properties and Case Histories. Bali, Indonesia, 27 - 48.\
\
Robertson, P. K. and Cabal, K. (2022). Guide to cone penetration testing. 7th Ed. Gregg Drilling LLC
\
\

## API Index

- link to [`find_cpt_column_names`](@ref)
- link to [`get_least_squares_interpolator`](@ref)
