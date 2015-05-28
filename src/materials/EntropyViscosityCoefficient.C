#include "EntropyViscosityCoefficient.h"

template<>
InputParameters validParams<EntropyViscosityCoefficient>()
{
  InputParameters params = validParams<Material>();

  // Coupled variables:
  params.addRequiredCoupledVar("radiation_flux", "Variable storing the radiation flux");
  // Jump values
  params.addCoupledVar("entropy_flux_jump", "variable storing the jump of the entropy flux across the surface");
  // Speed of light constant:
  params.addParam<Real>("c", 1., "speed of light value");
  // Angular:
  params.addParam<Real>("angular_direction", 1., "angular direction of the radiation: +1 or -1");
  // Cconstant parameter:
  params.addParam<Real>("Ce", 1., "Coefficient for entropy residual");
  params.addParam<Real>("Cj", 1., "Coefficient for jump");
  // Name of the pps computing the normalization parameter
  params.addRequiredParam<std::string>("name_pps_for_normalization", "Name of the pps computing the normalization parameter");

  return params;
}

EntropyViscosityCoefficient::EntropyViscosityCoefficient(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Coupled variables:
    _u_old(coupledValueOld("radiation_flux")),
    _u_older(coupledValueOlder("radiation_flux")),
    _grad_u_old(coupledGradientOld("radiation_flux")),
    // Jump values:
    _jump(isCoupled("entropy_flux_jump") ? coupledValue("entropy_flux_jump") : _zero),
    // Speed of light constant:
    _c(getParam<Real>("c")),
    // Angular
    _omega(getParam<Real>("angular_direction")),
    // Parameters
    _Ce(getParam<double>("Ce")),
    _Cj(getParam<double>("Cj")),
    // Name of the pps computing the normalization parameter
    _pps_name(getParam<std::string>("name_pps_for_normalization")),
    // Declare material properties storing entropy viscosity coefficient.
    _kappa(declareProperty<Real>("entropy_viscosity_coefficient")),
    // Get material: cross section
    _sigma(getMaterialProperty<Real>("sigma"))
{
  mooseAssert(_Cj<0, this->name() << ": the coefficient Cj has to be positive.");
}

void
EntropyViscosityCoefficient::computeQpProperties()
{
  // Cell size
  Real h_cell = std::pow(_current_elem->volume(),1./_mesh.dimension());

  // Compute the entropy values at time 'n' and 'n-1' and its derivative 'dEdu'
  Real E_old = 0.5*_u_old[_qp]*_u_old[_qp];
  Real E_older = 0.5*_u_older[_qp]*_u_older[_qp];
  Real dEdu_old = _u_old[_qp];

  // Compute the pressure contribution to the residual:
  Real residual = (E_old-E_older)/_dt + dEdu_old*(_omega*_grad_u_old[_qp](0)+_sigma[_qp]*_u_old[_qp]);

  // Get normalization parameter from pps
  Real norm = std::max(getPostprocessorValueByName(_pps_name), 1.e-6);

  // Return entropy viscosity coefficient:
  _kappa[_qp] = h_cell*h_cell*(_Ce*std::fabs(residual)+_Cj*_jump[_qp])/norm;
}