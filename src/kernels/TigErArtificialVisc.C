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

#include "TigErArtificialVisc.h"
/**
This function computes the artificial dissipative terms. It only works in 1-D
 */
template<>
InputParameters validParams<TigErArtificialVisc>()
{
  InputParameters params = validParams<Kernel>();

  return params;
}

TigErArtificialVisc::TigErArtificialVisc(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Material properties:
    _kappa(getMaterialProperty<Real>("kappa"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real TigErArtificialVisc::computeQpResidual()
{
  // Compute the viscosity term

  // Return the term
  return _kappa[_qp]*_grad_u[_qp]*_grad_test[_i][_qp];
}

Real TigErArtificialVisc::computeQpJacobian()
{
  return 0.;
}

Real TigErArtificialVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}