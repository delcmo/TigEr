#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST

###### Constants #######
c = 1.
angular_direction = 1.
[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0.
  xmax = 1.
  elem_type = EDGE2
[]

#############################################################################
#                             FUNCTIONS                                     #
#############################################################################

[Functions]
  [./ic_func]
    axis = 0
    type                      = PiecewiseLinear
    x                         = '0    0.5   0.50001  1.'
    y                         = '1. 1.  0.  0.'
#    x                         = '0  1.'
#    y                         = '2. 1.'
  [../]

  [./ic_func_sin]
    type                      = ParsedFunction
    value = sin(2*pi*x)
  [../]

  [./exact_sin]
    type                      = ParsedFunction
    value = sin(2*pi*x-t)
  [../]
[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################

[Variables]
  [./radiation]
    family = LAGRANGE
    order = FIRST
    scaling = 1e+0
      [./InitialCondition]
        type = FunctionIC
        function = ic_func
#        value = 10.
      [../]
   [../]
[]

############################################################################################################
#                                            KERNELS                                                       #
############################################################################################################

[Kernels]
  [./RadiationTime]
    type = TimeDerivative
    variable = radiation
    lumping = true
    implicit = true
  [../]

  [./RadiationAdvection]
    type = Diffusion # TigErAdvection
    variable = radiation
    implicit = false
  [../]

#  [./RadiationCrossSection]
#    type = TigErCrossSection
#    variable = radiation
#  [../]

# [./RadiationArtificial]
#    type = TigErArtificialVisc
#    variable = radiation
#    implicit = false
#  [../]
[]

##############################################################################################
#                                       AUXILARY VARIABLES                                   #
##############################################################################################

[AuxVariables]
[]

##############################################################################################
#                                       AUXILARY KERNELS                                     #
##############################################################################################

[AuxKernels]
[]

##############################################################################################
#                                       MATERIALS                                     #
##############################################################################################

[Materials]
  [./kappaMat]
    type = GenericConstantMaterial
    prop_names = 'kappa'
    block = 0
    prop_values = 1.
  [../]

  [./sigmaMat]
    type = GenericConstantMaterial
    prop_names = 'sigma'
    block = 0
    prop_values = 0.
  [../]
[]

##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[BCs]
  [./RadiationRight]
    type = DirichletBC
#    type = FunctionDirichletBC
    variable = radiation
    value = 1.
#    function = exact_sin
    boundary = 'left'
    implicit = false
  [../]
  
  [./RadiationLeft]
    type = DirichletBC
    variable = radiation
    value = 0.
    boundary = 'right'
    implicit = false
  [../]
[]
##############################################################################################
#                                  PRECONDITIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

#[Preconditioning]
#  active = ' '
#  [./FDP_Newton]
#   type = FDP
#    full = true
#    solve_type = 'NEWTON' # 'PJFNK'
#petsc_options = '-snes_mf_operator -snes_ksp_ew'
#petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
#petsc_options_value = '1.e-12       ds             ds'
#  [../]
#[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient
  scheme = 'explicit-euler'
  solve_type = 'LINEAR'

  start_time = 0.0
  num_steps = 20
  dt = 0.00005
#
#  type = Transient
#  scheme = 'bdf2'
#  scheme = 'explicit-euler'
#  solve_type = 'LINEAR'
#  end_time = 0.3
#  dt = 1.e-4
#  dtmin = 1e-9
#  l_tol = 1e-8
#  nl_rel_tol = 1e-10
#  nl_abs_tol = 1e-7
#  l_max_its = 10
#  nl_max_its = 10
#  num_steps = 1 # 5000000
#  [./TimeStepper]
#    type = FunctionDT
#    time_t =  '0.     1.'
#    time_dt = '1.e-4  1.e-4'
#  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Outputs]
  output_initial = true
  interval = 1
  exodus = true
[]

##############################################################################################
#                                        DEBUG                                               #
##############################################################################################

#[Debug]
#  show_var_residual_norms = true
#[]
