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
    type = ParsedFunction
    value = 0
  [../]

  [./forcing_fn]
    type = ParsedFunction
#    value = x
    value = x+t
  [../]

  [./exact_fn]
    type = ParsedFunction
    value = t*x
  [../]
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE

    [./InitialCondition]
      type = FunctionIC
      function = ic
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

  [./ffn]
    type = UserForcingFunction
    variable = u
    function = forcing_fn
    implicit = false
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
  active = 'all'

  [./all]
    type = FunctionDirichletBC
    variable = u
    boundary = '0 1'
    function = exact_fn
    implicit = true
  [../]
[]

[Postprocessors]
  [./l2_err]
    type = ElementL2Error
    variable = u
    function = exact_fn
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'explicit-euler'
  solve_type = 'LINEAR'

  start_time = 0.0
  num_steps = 20
  dt = 0.00005
[]

[Outputs]
  output_initial = true
  exodus = true
  print_perf_log = true
  [./console]
    type = Console
    max_rows = 10
  [../]
[]
