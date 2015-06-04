#include "TigerApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

// time integrators
#include "ExplicitEulerFCT.h"

// kernels
#include "LumpedTimeDerivative.h"
#include "TigErAdvection.h"
#include "TigErCrossSection.h"
#include "TigErArtificialVisc.h"
#include "TigErMassMatrixDiffusion.h"
#include "TigErAntiDiffusionTerm.h"

// auxkernels
#include "SttResidualAux.h"
#include "Entropy.h"

// materials
#include "EntropyViscosityCoefficient.h"

// userobjects
#include "GetExtremumValueFromNeighbors.h"

// postprocessors
#include "InfiniteNormFromAverageValue.h"

template<>
InputParameters validParams<TigerApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

TigerApp::TigerApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  TigerApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  TigerApp::associateSyntax(_syntax, _action_factory);
}

TigerApp::~TigerApp()
{
}

extern "C" void TigerApp__registerApps() { TigerApp::registerApps(); }
void
TigerApp::registerApps()
{
  registerApp(TigerApp);
}

void
TigerApp::registerObjects(Factory & factory)
{
  // time integrator
  registerTimeIntegrator(ExplicitEulerFCT);

  // kernels
  registerKernel(LumpedTimeDerivative);
  registerKernel(TigErAdvection);
  registerKernel(TigErCrossSection);
  registerKernel(TigErArtificialVisc);
  registerKernel(TigErMassMatrixDiffusion);
  registerKernel(TigErAntiDiffusionTerm);

  // auxkernels
  registerAux(SttResidualAux);
  registerAux(Entropy);

  // materials
  registerMaterial(EntropyViscosityCoefficient);

  // userobjects
  registerUserObject(GetExtremumValueFromNeighbors);

  // postprocessors
  registerPostprocessor(InfiniteNormFromAverageValue);
}

void
TigerApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
