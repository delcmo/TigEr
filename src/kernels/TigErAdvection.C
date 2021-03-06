/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                               */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "TigErAdvection.h"

/**
This Kernel computes the convection flux of the continuity equation.
*/
template<>
InputParameters validParams<TigErAdvection>()
{
  InputParameters params = validParams<Kernel>();

  // Speed of light constant:
  params.addParam<Real>("c", 1., "speed of light value");
  // Angular:
  params.addParam<Real>("angular_direction", 1., "angular direction of the radiation: +1 or -1");
  
  return params;
}

TigErAdvection::TigErAdvection(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Speed of light constant:
    _c(getParam<Real>("c")),
    // Angular
    _omega(getParam<Real>("angular_direction"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real TigErAdvection::computeQpResidual()
{
  // Return advection term
  return _c*_omega*_grad_u[_qp](0)*_test[_i][_qp];
//  return -_c*_omega*_u[_qp]*_grad_test[_i][_qp](0);
}

Real TigErAdvection::computeQpJacobian()
{
  return 0.;
}

Real TigErAdvection::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.;
}