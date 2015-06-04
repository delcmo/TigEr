[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0.
  xmax = 1
  nx = 10
[]

[Functions]
  [./ic1]
    axis = 0
    type                      = PiecewiseLinear
    x                         = '0    0.5   0.51  1.'
    y                         = '0. 0.  0.  0.'
  [../]
  
  [./ic2]
    axis = 0
    type                      = PiecewiseLinear
    x                         = '0  0.3  0.4  0.6   0.8  0.95 1.'
    y                         = '1. 3.   3.   2     2    0.   0.'
  [../]
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE

    [./InitialCondition]
      type = FunctionIC
      function = ic1
    [../]
  [../]
[]

[Kernels]
  [./ie]
    type = LumpedTimeDerivative
    variable = u
    implicit = true
  [../]

  [./advection]
    type = TigErAdvection
    variable = u
    implicit = false
  [../]

  [./diff]
    type = TigErArtificialVisc
    variable = u
    implicit = false
  [../]
  
  [./antidiffusionterm]
    type = TigErAntiDiffusionTerm
    variable = u
    max_nodal_values = u_max_node
    min_nodal_values = u_min_node
    implicit = false
  [../]
[]

[AuxVariables]
  [./entropy]
    order = FIRST
    family = LAGRANGE
  [../]

  [./u_max_node]
    order = FIRST
    family = LAGRANGE
  [../]
  
  [./u_min_node]
    order = FIRST
    family = LAGRANGE
  [../]
  
  [./kappa_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./EntropyAK]
    type = Entropy
    variable = entropy
    radiation_flux = u
  [../]

  [./kappaAK]
    type = MaterialRealAux
    variable = kappa_aux
    property = entropy_viscosity_coefficient
  [../]
[]

[UserObjects]
  [./u_max_uo]
    type = GetExtremumValueFromNeighbors
    variable = u
    variable_out = u_max_node
    execute_on = 'timestep_end'
  [../]
  
  [./u_min_uo]
    type = GetExtremumValueFromNeighbors
    variable = u
    variable_out = u_min_node
    compute_maximum = false
    execute_on = 'timestep_end'
  [../]
[]

[Postprocessors]
  [./average_entropy]
    type = ElementAverageValue
    variable = entropy
    execute_on = 'timestep_end'
  [../]
  
  [./infinite_norm]
    type = InfiniteNormFromAverageValue
    variable = entropy
    name_pps_average = average_entropy
    execute_on = 'timestep_end'
  [../]
[]

[Materials]
  [./sigmaMat]
    type = GenericConstantMaterial
    prop_names = 'sigma'
    block = 0
    prop_values = 0.
  [../]
  
#  [./kappaMat]
#    type = GenericConstantMaterial
#    prop_names = 'entropy_viscosity_coefficient'
#    block = 0
#    prop_values = 0.
#  [../]
  
  [./EntropyViscosityCoefficient]
    type = EntropyViscosityCoefficient
    radiation_flux = u
    name_pps_for_normalization = infinite_norm
    block = 0    
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = '0'
    value = 1.
    implicit = true
  [../]

 [./right]
   type = DirichletBC
   variable = u
   boundary = '1'
   value = 0.
   implicit = true
 [../]
[]

[Executioner]
  type = Transient
  scheme = 'explicit-euler-fct'
  solve_type = 'LINEAR'
  petsc_options = '-snes_converged_reason'

  start_time = 0.0
  end_time = 0.3
  num_steps = 1 # 5000
  dt = 0.0004

  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]
[]

[Outputs]
  output_initial = true
  exodus = true
  print_perf_log = true
  print_linear_residuals = true
  [./console]
    type = Console
  [../]
[]
