#ifndef ENTROPYVISCOSITYCOEFFICIENT_H
#define ENTROPYVISCOSITYCOEFFICIENT_H

#include "Material.h"
#include "MaterialProperty.h"

//Forward Declarations
class EntropyViscosityCoefficient;

template<>
InputParameters validParams<EntropyViscosityCoefficient>();

class EntropyViscosityCoefficient : public Material
{
public:
  EntropyViscosityCoefficient(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:

  // Coupled variables
  VariableValue & _u_old;
  VariableValue & _u_older;
  VariableGradient & _grad_u_old;

  // Jump value:
  VariableValue & _jump;

  // Speed of light constant
  Real _c;

  // Angular
  Real _omega;

  // Parameters
  Real _Ce;
  Real _Cj;

  // Name of the pps computing the normalization parameter
  std::string _pps_name;

  // Material property: entropy viscosity coefficient.
  MaterialProperty<Real> & _kappa;

  // Material property: cross section
  MaterialProperty<Real> & _sigma;
};

#endif // ENTROPYVISCOSITYCOEFFICIENT_H