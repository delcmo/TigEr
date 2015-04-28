/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
/**
This function computes the fluid internal energy 'rhoe' from the conservative variables. It is dimension agnostic.
**/
#include "SttResidualAux.h"

template<>
InputParameters validParams<SttResidualAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Speed of light value
  params.addParam<Real>("c", 1., "speed of light value");
  // Angular
  params.addParam<Real>("angular_direction", 1., "angular direction of the radiation: +1 or -1");
  // Coupled variables
  params.addRequiredCoupledVar("radiation", "variable that computes the radiation");
  
  return params;
}

SttResidualAux::SttResidualAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables:
    _radiation(coupledValue("radiation")),
    _rad_grad(coupledGradient("radiation")),
    // Cross section:
    _sigma(getMaterialProperty<Real>("sigma")),
    // Speed of light:
    _c(getParam<Real>("c")),
    // Angular
    _Omega(getParam<Real>("angular_direction"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real
SttResidualAux::computeValue()
{
  // Return stt value residual at quadrature point _qp:
  return _c*( _Omega*_rad_grad[_qp](0)+_sigma[_qp]*_radiation[_qp] );
}