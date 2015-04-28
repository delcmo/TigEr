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

#ifndef STTRESIDUALAUX_H
#define STTRESIDUALAUX_H

#include "AuxKernel.h"

class SttResidualAux;

template<>
InputParameters validParams<SttResidualAux>();

class SttResidualAux : public AuxKernel
{
public:

  SttResidualAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  // Coupled variables:
  VariableValue & _radiation;
  VariableGradient & _rad_grad;

  // Cross section (material property)
  MaterialProperty<Real> & _sigma;

  // Speed of light
  Real _c;

  // Angular
  Real _Omega;
};

#endif // STTRESIDUALAUX_H