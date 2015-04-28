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

#ifndef TIGERCROSSSECTION_H
#define TIGERCROSSSECTION_H

#include "Kernel.h"

// Forward Declarations
class TigErCrossSection;

template<>
InputParameters validParams<TigErCrossSection>();

class TigErCrossSection : public Kernel
{
public:

  TigErCrossSection(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:
  // Speed of light constant:
  Real _c;

  // Cross section (material property)
  MaterialProperty<Real> & _sigma;
};

#endif // TIGERCROSSSECTION_H
