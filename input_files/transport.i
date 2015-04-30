#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST
lumping = true
#implicit = false

###### Constants #######
c = 1.
angular_direction = 1.
[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 32
  xmin = 0.
  xmax = 1.
[]

#############################################################################
#                             FUNCTIONS                                     #
#############################################################################

[Functions]
  [./ic_func]
    axis = 0
    type                      = PiecewiseLinear
    x                         = '0    5   5.0001  10.'
    y                         = '1. 1.  0.  0.'
#    x                         = '0  1.'
#    y                         = '2. 1.'
  [../]
[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################

[Variables]
  [./radiation]
    family = LAGRANGE
    scaling = 1e+0
      [./InitialCondition]
        type = ConstantIC # FunctionIC
#        function = ic_func
        value = 0.
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
  [../]

  [./RadiationAdvection]
    type = TigErAdvection
    variable = radiation
  [../]

  [./RadiationCrossSection]
    type = TigErCrossSection
    variable = radiation
  [../]

  [./RadiationArtificial]
    type = TigErArtificialVisc
    variable = radiation
    avg_stt_res = low_order_visc
  [../]
[]

##############################################################################################
#                                       AUXILARY VARIABLES                                   #
##############################################################################################

[AuxVariables]
  [./low_order_visc]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

##############################################################################################
#                                       AUXILARY KERNELS                                     #
##############################################################################################

[AuxKernels]
  [./SttResidualAK]
    type = SttResidualAux
    variable = low_order_visc
    radiation = radiation
  [../]
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
    variable = radiation
    value = 1.
    boundary = 'left'
  [../]
  
  [./RadiationLeft]
    type = DirichletBC
    variable = radiation
    value = 0.
    boundary = 'right'
  [../]
[]
##############################################################################################
#                                  PRECONDITIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Preconditioning]
  active = 'FDP_Newton'
  [./FDP_Newton]
    type = FDP
    full = true
    solve_type = 'PJFNK'
#petsc_options = '-snes_mf_operator -snes_ksp_ew'
#petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
#petsc_options_value = '1.e-12       ds             ds'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient
  scheme = 'bdf2'
  end_time = 0.3
  dt = 1.e-4
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-7
  l_max_its = 10
  nl_max_its = 10
  num_steps = 5000000
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0.     1.'
    time_dt = '1.e-4  1.e-4'
  [../]
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
