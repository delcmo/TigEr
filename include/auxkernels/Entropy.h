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

#ifndef ENTROPY_H
#define ENTROPY_H

#include "AuxKernel.h"

class Entropy;

template<>
InputParameters validParams<Entropy>();

class Entropy : public AuxKernel
{
public:

  Entropy(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  // Coupled variables:
  VariableValue & _radiation;
};

#endif // ENTROPY_H