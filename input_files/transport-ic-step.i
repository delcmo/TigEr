[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0.
  xmax = 1
  nx = 200
  elem_type = EDGE2
[]

[Functions]
  [./ic]
    axis = 0
    type                      = PiecewiseLinear
    x                         = '0    0.5   0.51  1.'
    y                         = '1. 1.  0.  0.'
  [../]
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE

    [./InitialCondition]
      type = FunctionIC
      function = ic
      type = ConstantIC
      value = 0.
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
#    implicit = false
#   [./]
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

[Executioner]
  type = Transient
  scheme = 'explicit-euler'
  solve_type = 'LINEAR'
  petsc_options = '-snes_converged_reason'

  start_time = 0.0
  end_time = 0.3
#  num_steps = 1000 # 5000
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
