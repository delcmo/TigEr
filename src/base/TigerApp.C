#include "TigerApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

// kernels
#include "TigErAdvection.h"
#include "TigErCrossSection.h"
#include "TigErArtificialVisc.h"
#include "TigErMassMatrixDiffusion.h"

// auxkernels
#include "SttResidualAux.h"

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
  // kernels
  registerKernel(TigErAdvection);
  registerKernel(TigErCrossSection);
  registerKernel(TigErArtificialVisc);
  registerKernel(TigErMassMatrixDiffusion);

  // auxkernels
  registerAux(SttResidualAux);
}

void
TigerApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
