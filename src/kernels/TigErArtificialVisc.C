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

  // Speed of light constant:
  params.addParam<Real>("c", 1., "speed of light value");
  // Angular:
  params.addParam<Real>("angular_direction", 1., "angular direction of the radiation: +1 or -1");
  // Elemental value of the steady-state residual
  params.addRequiredCoupledVar("avg_stt_res", "elemental average value of the steady-state residual");

  return params;
}

TigErArtificialVisc::TigErArtificialVisc(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Averaged value of the steady-state residual
    _avg_stt_res(coupledValue("avg_stt_res")),
    // Speed of light constant:
    _c(getParam<Real>("c")),
    // Angular
    _omega(getParam<Real>("angular_direction")),
    // Material properties:
    _kappa(getMaterialProperty<Real>("kappa")),
    _sigma(getMaterialProperty<Real>("sigma"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real TigErArtificialVisc::computeQpResidual()
{
//  std::cout<<"&&&&&&&&&&&&&&"<<std::endl;
  // Get the nodal values and the shape functions
  VariableValue var_nodal = _var.nodalValue();
  VariablePhiValue phi = _var.phi();
  VariablePhiGradient phi_grad = _var.gradPhi();
//  std::cout<<var_nodal.size()<<std::endl;
//  std::cout<<phi.size()<<std::endl;
//  std::cout<<_phi.size()<<std::endl;
//  std::cout<<"nodal value="<<var_nodal[0]<<" and "<<var_nodal[1]<<std::endl;
//  std::cout<<"qp value="<<_u[0]<<" and "<<_u[1]<<std::endl;
//  std::cout<<"qp value="<<phi[0][0]*var_nodal[0]+phi[1][0]*var_nodal[1]<<" and "<<_u[1]<<std::endl;
//  _var.isNodal()?std::cout<<"true"<<std::endl : std::cout<<"false"<<std::endl;

  // Compute A_{i,j}, component (i,j) of the steady-state matrix
  Real Aij[phi.size()][phi.size()];
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    {
      Aij[ivar][jvar] = 0.;
      for (unsigned int qp=0; qp<_qrule->n_points(); qp++)
      {
        Aij[ivar][jvar] += _omega*phi_grad[jvar][qp](0)+_sigma[_qp]*phi[jvar][qp];
        Aij[ivar][jvar] *= _c*_test[ivar][qp]*_coord[_qp]*_JxW[_qp];
      }
    }

  // Compute the bilinear form b_k(\phi_j, \phi_i)
  Real b_k[phi.size()][phi.size()];
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
      b_k[ivar][jvar] = jvar == ivar ? _current_elem_volume : -_current_elem_volume;

//  for (unsigned int jvar=0; jvar<phi.size(); jvar++)
//    for (unsigned int ivar=0; ivar<phi.size(); ivar++)
//    {
//      std::cout<<"Aij["<<ivar<<"]["<<jvar<<"]="<<Aij[jvar][ivar]<<std::endl;
//      std::cout<<"b_k["<<ivar<<"]["<<jvar<<"]="<<b_k[jvar][ivar]<<std::endl;
//    }

  // Compute the viscosity term
  Real sum_b_k(0.), max_Aij(0.);
  for (unsigned int ivar=0; ivar<phi.size(); ivar++)
    for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    {
      if (ivar != jvar)
      {
        max_Aij = std::max(max_Aij, std::max(0., Aij[ivar][jvar]/(-b_k[ivar][jvar])));
        sum_b_k += b_k[ivar][jvar];
      }
    }
//  std::cout<<"max_Aij="<<max_Aij<<std::endl;
//  std::cout<<"sum_b_k="<<sum_b_k<<std::endl;
  Real visc = max_Aij;///(-sum_b_k);
  
  if (visc<0)
    mooseError("The low-order viscosity coefficient computed in '"<<this->name()<<"' is locally negative.");
//  std::cout<<"visc="<<visc<<std::endl;
//  std::cout<<"weight="<<_coord[_qp]<<std::endl;

  // Compute the artificial diffusion term
  Real diff_term_laplace = _grad_u[_qp](0)*_grad_test[_i][_qp](0);
  Real diff_term(0.);
//  std::cout<<"diff_term="<<diff_term<<std::endl;
  for (unsigned int jvar=0; jvar<phi.size(); jvar++)
    diff_term += var_nodal[jvar]*b_k[_i][jvar];
  diff_term *= visc/(_JxW[_qp]*_coord[_qp]);

  //  std::cout<<"current_elem_volume="<<_current_elem_volume<<std::endl;
//  std::cout<<"diff_term="<<diff_term<<std::endl;
//  std::cout<<"diff_term_la="<<diff_term_laplace<<std::endl;

  // Return value
  return diff_term;
}

Real TigErArtificialVisc::computeQpJacobian()
{
  return 0.;
}

Real TigErArtificialVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}