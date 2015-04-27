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

#ifndef TIGERADVECTION_H
#define TIGERADVECTION_H

#include "Kernel.h"

class TigErAdvection;

template<>
InputParameters validParams<TigErAdvection>();
class TigErAdvection : public Kernel
{
public:

  TigErAdvection(const std::string & name,
             InputParameters parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:

  // Angular
  Real _omega;
};

#endif // TIGERADVECTION_H