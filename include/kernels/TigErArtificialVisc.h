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

#ifndef TIGERARTIFICIALVISC_H
#define TIGERARTIFICIALVISC_H

#include "Kernel.h"

// Forward Declarations
class TigErArtificialVisc;

template<>
InputParameters validParams<TigErArtificialVisc>();

class TigErArtificialVisc : public Kernel
{
public:

  TigErArtificialVisc(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:
  // Material property:
  MaterialProperty<Real> & _kappa;
};

#endif // TIGERARTIFICIALVISC_H
