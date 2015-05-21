[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0.
  xmax = 1
  nx = 100
#  elem_type = EDGE3
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
      function = ic2
#      type = ConstantIC
#      value = 0.
    [../]
  [../]
[]

[Kernels]
  [./ie]
    type = TimeDerivative
    variable = u
    lumping = true
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

#   [./lump_mass_matrix]
#    type = TigErMassMatrixDiffusion
#    variable = u
#    implicit = true
#   [./]
[]

[AuxVariables]
  [./u_max_node]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
# active = ''
  [./u_max_uo]
    type = GetExtremumValueFromNeighbors
    variable = u
    variable_out = u_max_node
    execute_on = 'timestep_begin'
  [../]
[]

[Materials]
  [./sigmaMat]
    type = GenericConstantMaterial
    prop_names = 'sigma'
    block = 0
    prop_values = 0.
  [../]
[]

[BCs]
  active = 'left right'

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

#[Postprocessors]
#  [./u_max_uo]
#   type = GetExtremumValueFromNeighbors
#    variable = u
#    variable_out = u_max_node
#    execute_on = 'timestep_begin'
#  [../]
#[]

[Executioner]
  type = Transient
  scheme = 'rk-2'
  solve_type = 'LINEAR'
  petsc_options = '-snes_converged_reason'

  start_time = 0.0
  end_time = 0.3
  num_steps = 10 # 5000
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
