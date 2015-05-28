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
#include "Entropy.h"

template<>
InputParameters validParams<Entropy>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variable
  params.addRequiredCoupledVar("radiation_flux", "variable that computes the radiation");
  
  return params;
}

Entropy::Entropy(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables:
    _radiation(coupledValue("radiation_flux"))
{
  mooseAssert(_mesh.dimension()!=1, "The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real
Entropy::computeValue()
{
  // Return the value of the entropy:
  return 0.5*_radiation[_qp]*_radiation[_qp];
}