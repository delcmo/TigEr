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

#include "TigErCrossSection.h"
/**
This function computes the cross section term for the 1-d transport equation. It only works in 1-D
 */
template<>
InputParameters validParams<TigErCrossSection>()
{
  InputParameters params = validParams<Kernel>();

  // Speed of light constant:
  params.addParam<Real>("c", 1., "speed of light value");

  return params;
}

TigErCrossSection::TigErCrossSection(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Speed of light constant:
    _c(getParam<Real>("c")),
    // Cross section:
    _sigma(getMaterialProperty<Real>("sigma"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real TigErCrossSection::computeQpResidual()
{
  // Return value:
  return _c*_sigma[_qp]*_u[_qp]*_test[_i][_qp];
}

Real TigErCrossSection::computeQpJacobian()
{
  return 0.;
}

Real TigErCrossSection::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}