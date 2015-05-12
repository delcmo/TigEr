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

#include "TigErMassMatrixDiffusion.h"
/**
This function computes the artificial dissipative terms. It only works in 1-D
 */
template<>
InputParameters validParams<TigErMassMatrixDiffusion>()
{
  InputParameters params = validParams<Kernel>();

  // Speed of light constant:
  params.addParam<Real>("c", 1., "speed of light value");

  return params;
}

TigErMassMatrixDiffusion::TigErMassMatrixDiffusion(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // get the nodal values
    _u_nodal(_is_implicit ? _var.nodalValue() : _var.nodalValueOld()),
    // Speed of light constant:
    _c(getParam<Real>("c"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real TigErMassMatrixDiffusion::computeQpResidual()
{
  // Return value
  return -_c*(_u[_qp]-_u_nodal[_i])*_test[_i][_qp];
}

Real TigErMassMatrixDiffusion::computeQpJacobian()
{
  return 0.;
}

Real TigErMassMatrixDiffusion::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}