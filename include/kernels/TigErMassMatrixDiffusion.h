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

#ifndef TIGERMASSMATRIXDIFFUSION_H
#define TIGERMASSMATRIXDIFFUSION_H

#include "Kernel.h"

// Forward Declarations
class TigErMassMatrixDiffusion;

template<>
InputParameters validParams<TigErMassMatrixDiffusion>();

class TigErMassMatrixDiffusion : public Kernel
{
public:

  TigErMassMatrixDiffusion(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:

  // Nodal values
  VariableValue & _u_nodal;

  // Speed of light constant
  Real _c;
};

#endif // TIGERMASSMATRIXDIFFUSION_H