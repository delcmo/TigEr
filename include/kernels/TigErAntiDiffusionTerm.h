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

#ifndef TIGERANTIDIFFUSIONTERM_H
#define TIGERANTIDIFFUSIONTERM_H

#include "Kernel.h"

// Forward Declarations
class TigErAntiDiffusionTerm;

template<>
InputParameters validParams<TigErAntiDiffusionTerm>();

class TigErAntiDiffusionTerm : public Kernel
{
public:

  TigErAntiDiffusionTerm(const std::string & name,
             InputParameters parameters);

protected:

  virtual void computeResidual();

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:

  // Nodal values
  VariableValue & _u_nodal_old;
  VariableValue & _u_nodal;

  // Coupled aux variables
  VariableValue & _U_plus;
  VariableValue & _U_minus;

  // Constants
  Real _c;
  Real _omega;

  // Material property:
  const MaterialProperty<Real> & _kappa;
  const MaterialProperty<Real> & _sigma;
};

#endif // TIGERANTIDIFFUSIONTERM_H