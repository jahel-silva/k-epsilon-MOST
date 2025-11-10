# k-epsilon-MOST
User-defined functions used in the dissertation entitled: Turbulence modeling the horizontally homogeneous atmospheric surface layer for wind energy applications using the k-epsilon model
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
UDF_K_E_FP_cstB_wake_CBLv5 - User defined function used in the wake analysis under stable condition. Includes the modifications needed to achieve horizontal homogeneity following Baungaard et al 2022, and the source term for the actuator disk.

UDF_K_E_FP_cstB_wake_SBLv5 - User defined function used in the wake analysis under unstable condition. Includes the modifications needed to achieve horizontal homogeneity following Baungaard et al 2022, and the source term for the actuator disk.

***IMPORTANT NOTES:
- UDFs were run in Ansys Fluent 19R2.
- The UDFs provided need to be compiled in the Ansys Fluent environment.
- The correct quantity of user-defined memories need to be set in Ansys environment. 
  Important to set the advanced features in Ansys Fluent:
  Type the following command in the console: 
  \solve\set\expert
  Enable the option "Keep temporary solver memory from being freed"
- The model includes, a modified wall function, a modified turbulent viscosity, a define_adjust macro, a define_initi macro, and the specified source terms. It is needed to hook the adjust (UDM_storage) and initialization functions (my_init_function). 
- The actuator disk source term is specific to the mesh developed and needs to be adapted to different mesh types.
