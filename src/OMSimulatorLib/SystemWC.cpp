/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
 * THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the Open Source Modelica
 * Consortium (OSMC) Public License (OSMC-PL) are obtained
 * from OSMC, either from the above address,
 * from the URLs: http://www.ida.liu.se/projects/OpenModelica or
 * http://www.openmodelica.org, and in the OpenModelica distribution.
 * GNU version 3 is obtained from: http://www.gnu.org/copyleft/gpl.html.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of  MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

#include "SystemWC.h"

#include "Component.h"
#include "ComponentFMUCS.h"
#include "Flags.h"
#include "Model.h"
#include "ssd/Tags.h"
#include "SystemTLM.h"
#include "Types.h"
#include "math.h"
oms3::SystemWC::SystemWC(const ComRef& cref, Model* parentModel, System* parentSystem)
  : oms3::System(cref, oms_system_wc, parentModel, parentSystem)
{
}

oms3::SystemWC::~SystemWC()
{
  if (derBuffer)
    delete[] derBuffer;
}

oms3::System* oms3::SystemWC::NewSystem(const oms3::ComRef& cref, oms3::Model* parentModel, oms3::System* parentSystem)
{
  if (!cref.isValidIdent())
  {
    logError_InvalidIdent(cref);
    return NULL;
  }

  if ((parentModel && parentSystem) || (!parentModel && !parentSystem))
  {
    logError_InternalError;
    return NULL;
  }

  System* system = new SystemWC(cref, parentModel, parentSystem);
  return system;
}

oms_status_enu_t oms3::SystemWC::exportToSSD_SimulationInformation(pugi::xml_node& node) const
{
  pugi::xml_node node_simulation_information = node.append_child(oms2::ssd::ssd_simulation_information);

  pugi::xml_node node_solver = node_simulation_information.append_child("FixedStepMaster");
  node_solver.append_attribute("description") = solverName.c_str();
  node_solver.append_attribute("stepSize") = std::to_string(stepSize).c_str();

  return oms_status_ok;
}

oms_status_enu_t oms3::SystemWC::importFromSSD_SimulationInformation(const pugi::xml_node& node)
{
  solverName = node.child("FixedStepMaster").attribute("description").as_string();
  stepSize = node.child("FixedStepMaster").attribute("stepSize").as_double();
  return oms_status_ok;
}

oms_status_enu_t oms3::SystemWC::instantiate()
{
  time = getModel()->getStartTime();

  for (const auto& subsystem : getSubSystems())
    if (oms_status_ok != subsystem.second->instantiate())
      return oms_status_error;

  for (const auto& component : getComponents())
    if (oms_status_ok != component.second->instantiate())
      return oms_status_error;

  return oms_status_ok;
}

unsigned int oms3::SystemWC::getMaxOutputDerivativeOrder()
{
  unsigned int order = 0;

  for (auto& component : getComponents())
  {
    if (oms_component_fmu == component.second->getType())
      if (order < component.second->getFMUInfo()->getMaxOutputDerivativeOrder())
        order = component.second->getFMUInfo()->getMaxOutputDerivativeOrder();
  }

  return order;
}

oms_status_enu_t oms3::SystemWC::initialize()
{
  clock.reset();
  CallClock callClock(clock);

  if (oms_status_ok != updateDependencyGraphs())
    return oms_status_error;

  if (oms_status_ok != updateInputs(initialUnknownsGraph))
    return oms_status_error;

  for (const auto& subsystem : getSubSystems())
    if (oms_status_ok != subsystem.second->initialize())
      return oms_status_error;

  for (const auto& component : getComponents())
    if (oms_status_ok != component.second->initialize())
      return oms_status_error;

  if (derBuffer)
    delete[] derBuffer;
  derBuffer = NULL;
  if (Flags::InputDerivatives())
    derBuffer = new double[getMaxOutputDerivativeOrder()];

  return oms_status_ok;
}

oms_status_enu_t oms3::SystemWC::terminate()
{
  for (const auto& subsystem : getSubSystems())
    if (oms_status_ok != subsystem.second->terminate())
      return oms_status_error;

  for (const auto& component : getComponents())
    if (oms_status_ok != component.second->terminate())
      return oms_status_error;

  return oms_status_ok;
}

oms_status_enu_t oms3::SystemWC::reset()
{
  for (const auto& subsystem : getSubSystems())
    if (oms_status_ok != subsystem.second->reset())
      return oms_status_error;

  for (const auto& component : getComponents())
    if (oms_status_ok != component.second->reset())
      return oms_status_error;

  return oms_status_ok;
}

oms_status_enu_t oms3::SystemWC::stepUntil(double stopTime, void (*cb)(const char* ident, double time, oms_status_enu_t status))
{
  CallClock callClock(clock);
  ComRef modelName = this->getModel()->getCref();

  fmi2_status_t fmi_status;
  double startTime=time;

  if (Flags::ProgressBar())
    logInfo("stepUntil [" + std::to_string(startTime) + "; " + std::to_string(stopTime) + "]");

  if(oms3::Flags::VariableStep()) // If variable step.
  {
    int fmuIndex = 0;
    std::map<ComRef, Component*> FMUcomponents;
    std::map<ComRef, Component*> canGetAndSetStateFMUcomponents;
    std::map<ComRef, Component*> noneFMUcomponents;
    bool firstTime = true;
    while (time < stopTime)
    {
      double tNext = time+stepSize;
      if (tNext > stopTime)
      {
        tNext = stopTime;
        stepSize = tNext-time;
      }
      logDebug("DEBUGGING: doStep: " + std::to_string(time) + " -> " + std::to_string(tNext));

      if(firstTime) // first time get a list of all components set up for which can get/Set fmu states.
      {
        firstTime = false;
        for (const auto& component : getComponents())
        {
          if (oms_component_fmu == component.second->getType()) // Check that its an FMU
          {
            if(dynamic_cast<ComponentFMUCS*>(component.second)->getFMUInfo()->getCanGetAndSetFMUstate())
            {
              canGetAndSetStateFMUcomponents.insert(std::pair<ComRef, Component*>(component.first,component.second));
            }
            else
            {
              FMUcomponents.insert(std::pair<ComRef, Component*>(component.first,component.second));
            }
          }
          else
          {
            noneFMUcomponents.insert(std::pair<ComRef, Component*>(component.first,component.second));
          }
        }
        if(canGetAndSetStateFMUcomponents.size() == 0) return oms_status_error;
      }
      oms_status_enu_t status;
      for (const auto& component : canGetAndSetStateFMUcomponents) // Get states and stepUntil for FMUs that can get state.
      {
        fmi2_import_t* fmu_in;
        fmi2_FMU_state_t s = NULL;
        fmu_in = dynamic_cast<ComponentFMUCS*>(component.second)->getFMU();
        s = NULL;
        fmi_status = fmi2_import_get_fmu_state(fmu_in,&s);
        sVect.push_back(s);
        fmiImportVect.push_back(fmu_in);

        status = component.second->stepUntil(tNext);
        if (oms_status_ok != status)
        {
          if (cb)
            cb(modelName.c_str(), tNext, status);
          return status;
        }
      }
      for (const auto& component : noneFMUcomponents) // stepUntil for noneFmus
      {
        status = component.second->stepUntil(tNext);
        if (oms_status_ok != status)
        {
          if (cb)
            cb(modelName.c_str(), tNext, status);
          return status;
        }
      }

      for (const auto& subsystem : getSubSystems()) // stepUntill for subsystems of ME FMUs, TODO: Fix rollback here too.
      {
        status = subsystem.second->stepUntil(tNext, NULL);
        if (oms_status_ok != status)
        {
          if (cb)
            cb(modelName.c_str(), tNext, status);
          return status;
        }
      }

      // get inputs at the end of all steps.
      std::vector<double> inputVect;
      std::vector<double> outputVect;

      if(oms_status_ok != getInput(outputsGraph,inputVect,FMUcomponents)) return oms_status_error;
      time = tNext;
      if (isTopLevelSystem())
        getModel()->emit(time);
      if(!FMUcomponents.empty()) // If we have FMUs that cant handle get/Set fmu state, we need to do the simulation order a bit different.
      {
        updateCanGetFMUs(outputsGraph,canGetAndSetStateFMUcomponents);
      }
      if(FMUcomponents.empty()) //If all FMUS can get and set states, we can do the error control after algebraic loops! YEY
      {        
        updateInputs(outputsGraph);
      }
      if (isTopLevelSystem())
        getModel()->emit(time);
      // get outputs after communication steps
      if(oms_status_ok != getOutput(outputsGraph,outputVect,FMUcomponents)) return oms_status_error;
        
      if (inputVect.size() != outputVect.size()) return oms_status_error;
      double biggestDifferance = 0.0;
      for (int n=0; n < inputVect.size();n++) // Calculate error in the FMUs we do error_control on.
      {
        double error;
        error = inputVect[n]-outputVect[n];
        if (error < 0) error = -error;
        if (error > biggestDifferance) biggestDifferance = error;
        logDebug("DEBUGGING: Input is: " + std::to_string(inputVect[n]) + " output is: " + std::to_string(outputVect[n]));
        logDebug("DEBUGGING: error is: " + std::to_string(error) + " Biggest Differance is: " + std::to_string(biggestDifferance));
      }
      double changeTolerance = 1e-2;
      double safety_factor = 0.90;
      double fixRatio = safety_factor*changeTolerance/biggestDifferance;
      logDebug("DEBUGGING: fixRatio is: " + std::to_string(fixRatio));

      if (biggestDifferance > changeTolerance) //Going to rollback.
      {
        logDebug("DEBUGGING: time is: " + std::to_string(time));
        logDebug("DEBUGGING: h is: " + std::to_string(stepSize));
        //Fix fmus
        for (int i=0; i<fmiImportVect.size(); ++i) // Reset all FMU states
        {
          fmi_status = fmi2_import_set_fmu_state(fmiImportVect[i], sVect[i]);
        }
        //Fix time
        logDebug("DEBUGGING: tNext is: " + std::to_string(tNext));
        time = tNext-stepSize;
        logDebug("DEBUGGING: time is after change: " + std::to_string(time));
        for (const auto& component : getComponents())
        {
          if (oms_component_fmu == component.second->getType())
          {
            dynamic_cast<ComponentFMUCS*>(component.second)->setFmuTime(time);
          }
        }
        //Fix Steptime.
        stepSize = stepSize*pow(fixRatio,0.25);
        logDebug("DEBUGGING: rollbacking, new h: " + std::to_string(stepSize));
      }
      else // Not going to rollback.
      {
        if (safety_factor*changeTolerance/biggestDifferance > 1)
        {
          stepSize = stepSize*pow(fixRatio,0.25);
          logDebug("DEBUGGING: Not rollbacking, new h: " + std::to_string(stepSize));
        }
        else
          logDebug("DEBUGGING: Not rollbacking, keeping old h: " + std::to_string(stepSize));
        if (!FMUcomponents.empty())
        {
          for (const auto& component : FMUcomponents) // These FMUs cant rollback, so only simulating them when we can.
          {
            status = component.second->stepUntil(tNext);
            if (oms_status_ok != status)
            {
              if (cb)
                cb(modelName.c_str(), tNext, status);
              return status;
            }           
          }    
          if (isTopLevelSystem())
            getModel()->emit(time); 
          updateCantGetFMUs(outputsGraph,FMUcomponents); // Update the last FMUS
          if (isTopLevelSystem())
            getModel()->emit(time);        
        }
      }
           
      fmiImportVect.clear();
      sVect.clear();
      if (cb)
        cb(modelName.c_str(), time, oms_status_ok);

      if (Flags::ProgressBar())
        Log::ProgressBar(startTime, stopTime, time);

      if (isTopLevelSystem() && getModel()->cancelSimulation())
      {
        cb(modelName.c_str(), time, oms_status_discard);
        return oms_status_discard;
      }
      
    }
  }
  else // If fixed step.
  {
    while (time < stopTime)
    {
      double tNext = time+stepSize;
      if (tNext > stopTime)
        tNext = stopTime;

      logDebug("doStep: " + std::to_string(time) + " -> " + std::to_string(tNext));

      oms_status_enu_t status;
      for (const auto& subsystem : getSubSystems())
      {
        status = subsystem.second->stepUntil(tNext, NULL);
        if (oms_status_ok != status)
        {
          if (cb)
            cb(modelName.c_str(), tNext, status);
          return status;
        }
      }

      for (const auto& component : getComponents())
      {
        status = component.second->stepUntil(tNext);
        if (oms_status_ok != status)
        {
          if (cb)
            cb(modelName.c_str(), tNext, status);
          return status;
        }
      }

      time = tNext;
      if (isTopLevelSystem())
        getModel()->emit(time);
      updateInputs(outputsGraph);
      if (isTopLevelSystem())
        getModel()->emit(time);

      if (cb)
        cb(modelName.c_str(), time, oms_status_ok);

      if (Flags::ProgressBar())
        Log::ProgressBar(startTime, stopTime, time);

      if (isTopLevelSystem() && getModel()->cancelSimulation())
      {
        cb(modelName.c_str(), time, oms_status_discard);
        return oms_status_discard;
      }
    }
  }
  if (Flags::ProgressBar())
  {
    Log::ProgressBar(startTime, stopTime, time);
    Log::TerminateBar();
  }
  return oms_status_ok;
}

oms_status_enu_t oms3::SystemWC::getRealOutputDerivative(const ComRef& cref, double*& value, unsigned int& order)
{
  if (!value)
    return oms_status_ok;

  switch (getModel()->getModelState())
  {
    case oms_modelState_initialization:
      return oms_status_error;
    case oms_modelState_instantiated:
    case oms_modelState_simulation:
      break;
    default:
      return logError_ModelInWrongState(getModel());
  }

  oms3::ComRef tail(cref);
  oms3::ComRef head = tail.pop_front();

  auto component = getComponents().find(head);
  if (component != getComponents().end() && oms_component_fmu == component->second->getType())
  {
    order = component->second->getFMUInfo()->getMaxOutputDerivativeOrder();
    if (order > 0)
      return dynamic_cast<ComponentFMUCS*>(component->second)->getRealOutputDerivative(tail, value);
  }

  return oms_status_error;
}

oms_status_enu_t oms3::SystemWC::setRealInputDerivative(const ComRef& cref, double* value, unsigned int order)
{
  if (!value)
    return oms_status_ok;

  switch (getModel()->getModelState())
  {
    case oms_modelState_instantiated:
    case oms_modelState_initialization:
    case oms_modelState_simulation:
      break;
    default:
      return logError_ModelInWrongState(getModel());
  }

  oms3::ComRef tail(cref);
  oms3::ComRef head = tail.pop_front();

  auto component = getComponents().find(head);
  if (component != getComponents().end() && oms_component_fmu == component->second->getType())
  {
    if (order > 0)
      return dynamic_cast<ComponentFMUCS*>(component->second)->setRealInputDerivative(tail, value, order);
  }

  return oms_status_error;
}



oms_status_enu_t oms3::SystemWC::getInput(oms3::DirectedGraph& graph, std::vector<double>& inputVect,std::map<ComRef, Component*> FMUcomponents)
{
  const std::vector< std::vector< std::pair<int, int> > >& sortedConnections = graph.getSortedConnections();
  inputVect.clear();
  if (!FMUcomponents.empty())
  {
    for(int i=0; i<sortedConnections.size(); i++)
    {
      int input = sortedConnections[i][0].second;
      oms3::ComRef inputName(graph.getNodes()[input].getName());
      oms3::ComRef inputModel = inputName.pop_front();
      int output = sortedConnections[i][0].first;
      oms3::ComRef outputName(graph.getNodes()[output].getName());
      oms3::ComRef outputModel = outputName.pop_front();
      if (FMUcomponents.find(inputModel)->first == inputModel || FMUcomponents.find(outputModel)->first == outputModel)
      {
        if (graph.getNodes()[input].getType() == oms_signal_type_real)
        {
          double inValue = 0.0;
          if (oms_status_ok != getReal(graph.getNodes()[input].getName(), inValue)) return oms_status_error;
          inputVect.push_back(inValue);
        }
      }
    }
  }
  else
  {
    for(int i=0; i<sortedConnections.size(); i++)
    {
      int input = sortedConnections[i][0].second;
      if (graph.getNodes()[input].getType() == oms_signal_type_real)
      {
        double inValue = 0.0;
        if (oms_status_ok != getReal(graph.getNodes()[input].getName(), inValue)) return oms_status_error;
        inputVect.push_back(inValue);
      }
    }
  }
  return oms_status_ok;
}

oms_status_enu_t oms3::SystemWC::getOutput(oms3::DirectedGraph& graph, std::vector<double>& outputVect,std::map<ComRef, Component*> FMUcomponents)
{
  const std::vector< std::vector< std::pair<int, int> > >& sortedConnections = graph.getSortedConnections();
  outputVect.clear();
  if (!FMUcomponents.empty())
  {
    for(int i=0; i<sortedConnections.size(); i++)
    {
      int output = sortedConnections[i][0].first;
      int input = sortedConnections[i][0].second;
      oms3::ComRef outputName(graph.getNodes()[output].getName());
      oms3::ComRef inputName(graph.getNodes()[input].getName());
      oms3::ComRef outputModel = outputName.pop_front();
      oms3::ComRef inputModel = inputName.pop_front();
      if (FMUcomponents.find(inputModel)->first == inputModel || FMUcomponents.find(outputModel)->first == outputModel)
      {
        if (graph.getNodes()[output].getType() == oms_signal_type_real)
        {
          double outValue = 0.0;
          if (oms_status_ok != getReal(graph.getNodes()[output].getName(), outValue)) return oms_status_error;
          outputVect.push_back(outValue);
        }
      }
    }
  }
  else
  {
    for(int i=0; i<sortedConnections.size(); i++)
    {
      int output = sortedConnections[i][0].first;
      if (graph.getNodes()[output].getType() == oms_signal_type_real)
      {
        double outValue = 0.0;
        if (oms_status_ok != getReal(graph.getNodes()[output].getName(), outValue)) return oms_status_error;
        outputVect.push_back(outValue);
      }
    }
  }
  return oms_status_ok;
}

oms_status_enu_t oms3::SystemWC::updateCanGetFMUs(oms3::DirectedGraph& graph,std::map<ComRef, Component*> FMUcomponents)
{
  const std::vector< std::vector< std::pair<int, int> > >& sortedConnections = graph.getSortedConnections();
  for(int i=0; i<sortedConnections.size(); i++)
  {
    if (sortedConnections[i].size() == 1)
    {
      int output = sortedConnections[i][0].first;
      int input = sortedConnections[i][0].second;
      oms3::ComRef outputName(graph.getNodes()[output].getName());
      oms3::ComRef inputName(graph.getNodes()[input].getName());
      oms3::ComRef outputModel = outputName.pop_front();
      oms3::ComRef inputModel = inputName.pop_front();
      if (FMUcomponents.find(inputModel)->first == inputModel && FMUcomponents.find(outputModel)->first == outputModel)
      {
        if (graph.getNodes()[input].getType() == oms_signal_type_real)
        {
          double value = 0.0;
          logDebug("DEBUGGING: We found a link between two fmus that can handle get and set");
          if (oms_status_ok != getReal(graph.getNodes()[output].getName(), value)) return oms_status_error;
          if (oms_status_ok != setReal(graph.getNodes()[input].getName(), value)) return oms_status_error;

          // derivatives
          if (derBuffer)
          {
            unsigned int order;
            if (oms_status_ok == getRealOutputDerivative(graph.getNodes()[output].getName(), derBuffer, order))
            {
              //logInfo(graph.getNodes()[output].getName() + " -> " + graph.getNodes()[input].getName() + ": " + std::to_string(derBuffer[0]));
              if (oms_status_ok != setRealInputDerivative(graph.getNodes()[input].getName(), derBuffer, order)) return oms_status_error;
            }
          }
        }
        else if (graph.getNodes()[input].getType() == oms_signal_type_integer)
        {
          int value = 0.0;
          if (oms_status_ok != getInteger(graph.getNodes()[output].getName(), value)) return oms_status_error;
          if (oms_status_ok != setInteger(graph.getNodes()[input].getName(), value)) return oms_status_error;
        }
        else if (graph.getNodes()[input].getType() == oms_signal_type_boolean)
        {
          bool value = 0.0;
          if (oms_status_ok != getBoolean(graph.getNodes()[output].getName(), value)) return oms_status_error;
          if (oms_status_ok != setBoolean(graph.getNodes()[input].getName(), value)) return oms_status_error;
        }
        else
          return logError_InternalError;

      }
    }
    else
    {
      bool OK = true;
      for (int j = 0; j < sortedConnections[j].size(); i++)
      {
        int output = sortedConnections[i][j].first;
        int input = sortedConnections[i][j].second;
        oms3::ComRef outputName(graph.getNodes()[output].getName());
        oms3::ComRef inputName(graph.getNodes()[input].getName());
        oms3::ComRef outputModel = outputName.pop_front();
        oms3::ComRef inputModel = inputName.pop_front();
        if (FMUcomponents.find(inputModel)->first == inputModel && FMUcomponents.find(outputModel)->first == outputModel)
        {
          OK = OK;
        }
        else
        {
          OK = false;
        }
      }
      if(OK)
      {
        if (oms_status_ok != solveAlgLoop(graph, sortedConnections[i])) return oms_status_error;  // This loop will only connect FMUs that can get state, and thus we can solve it.
      }
      else
      {
        logError("Algebraic Loop connecting FMUs that can get states and FMUs that can't get states, this will not work out. Restart the simulation with fixed step solver.");
        return oms_status_error;
      }
    }
  }
  return oms_status_ok;
}
oms_status_enu_t oms3::SystemWC::updateCantGetFMUs(oms3::DirectedGraph& graph,std::map<ComRef, Component*> FMUcomponents)
{
  const std::vector< std::vector< std::pair<int, int> > >& sortedConnections = graph.getSortedConnections();
  for(int i=0; i<sortedConnections.size(); i++)
  {
    if (sortedConnections[i].size() == 1)
    {
      int output = sortedConnections[i][0].first;
      int input = sortedConnections[i][0].second;
      oms3::ComRef outputName(graph.getNodes()[output].getName());
      oms3::ComRef inputName(graph.getNodes()[input].getName());
      oms3::ComRef outputModel = outputName.pop_front();
      oms3::ComRef inputModel = inputName.pop_front();
      if (FMUcomponents.find(inputModel)->first == inputModel || FMUcomponents.find(outputModel)->first == outputModel)
      {
        if (graph.getNodes()[input].getType() == oms_signal_type_real)
        {
          double value = 0.0;
          logDebug("DEBUGGING: We found a link between two fmus where atleast one cant handle get and set, so lets update this");
          if (oms_status_ok != getReal(graph.getNodes()[output].getName(), value)) return oms_status_error;
          if (oms_status_ok != setReal(graph.getNodes()[input].getName(), value)) return oms_status_error;

          // derivatives
          if (derBuffer)
          {
            unsigned int order;
            if (oms_status_ok == getRealOutputDerivative(graph.getNodes()[output].getName(), derBuffer, order))
            {
              //logInfo(graph.getNodes()[output].getName() + " -> " + graph.getNodes()[input].getName() + ": " + std::to_string(derBuffer[0]));
              if (oms_status_ok != setRealInputDerivative(graph.getNodes()[input].getName(), derBuffer, order)) return oms_status_error;
            }
          }
        }

        else if (graph.getNodes()[input].getType() == oms_signal_type_integer)
        {
          int value = 0.0;
          if (oms_status_ok != getInteger(graph.getNodes()[output].getName(), value)) return oms_status_error;
          if (oms_status_ok != setInteger(graph.getNodes()[input].getName(), value)) return oms_status_error;
        }
        else if (graph.getNodes()[input].getType() == oms_signal_type_boolean)
        {
          bool value = 0.0;
          if (oms_status_ok != getBoolean(graph.getNodes()[output].getName(), value)) return oms_status_error;
          if (oms_status_ok != setBoolean(graph.getNodes()[input].getName(), value)) return oms_status_error;
        }
        else
          return logError_InternalError;
      }

    }
    else
    {
      if (oms_status_ok != solveAlgLoop(graph, sortedConnections[i])) return oms_status_error; // This loop will only connect FMUs that cant get state, and thus we can solve it.
    }
  }
  return oms_status_ok;
}

oms_status_enu_t oms3::SystemWC::updateInputs(oms3::DirectedGraph& graph)
{
  CallClock callClock(clock);

  // input := output
  const std::vector< std::vector< std::pair<int, int> > >& sortedConnections = graph.getSortedConnections();
  for(int i=0; i<sortedConnections.size(); i++)
  {
    if (sortedConnections[i].size() == 1)
    {
      int output = sortedConnections[i][0].first;
      int input = sortedConnections[i][0].second;

      if (graph.getNodes()[input].getType() == oms_signal_type_real)
      {
        double value = 0.0;
        if (oms_status_ok != getReal(graph.getNodes()[output].getName(), value)) return oms_status_error;
        if (oms_status_ok != setReal(graph.getNodes()[input].getName(), value)) return oms_status_error;

        // derivatives
        if (derBuffer)
        {
          unsigned int order;
          if (oms_status_ok == getRealOutputDerivative(graph.getNodes()[output].getName(), derBuffer, order))
          {
            //logInfo(graph.getNodes()[output].getName() + " -> " + graph.getNodes()[input].getName() + ": " + std::to_string(derBuffer[0]));
            if (oms_status_ok != setRealInputDerivative(graph.getNodes()[input].getName(), derBuffer, order)) return oms_status_error;
          }
        }
      }
      else if (graph.getNodes()[input].getType() == oms_signal_type_integer)
      {
        int value = 0.0;
        if (oms_status_ok != getInteger(graph.getNodes()[output].getName(), value)) return oms_status_error;
        if (oms_status_ok != setInteger(graph.getNodes()[input].getName(), value)) return oms_status_error;
      }
      else if (graph.getNodes()[input].getType() == oms_signal_type_boolean)
      {
        bool value = 0.0;
        if (oms_status_ok != getBoolean(graph.getNodes()[output].getName(), value)) return oms_status_error;
        if (oms_status_ok != setBoolean(graph.getNodes()[input].getName(), value)) return oms_status_error;
      }
      else
        return logError_InternalError;
    }
    else
    {
      if (oms_status_ok != solveAlgLoop(graph, sortedConnections[i])) return oms_status_error;
    }
  }
  return oms_status_ok;
}


oms_status_enu_t oms3::SystemWC::solveAlgLoop(DirectedGraph& graph, const std::vector< std::pair<int, int> >& SCC)
{
  CallClock callClock(clock);

  const int size = SCC.size();
  const int maxIterations = 10;
  double maxRes;
  double *res = new double[size]();

  int it=0;
  do
  {
    it++;
    // get old values
    for (int i=0; i<size; ++i)
    {
      int output = SCC[i].first;
      if (oms_status_ok != getReal(graph.getNodes()[output].getName(), res[i]))
      {
        delete[] res;
        return oms_status_error;
      }
    }

    // update inputs
    for (int i=0; i<size; ++i)
    {
      int input = SCC[i].second;
      if (oms_status_ok != setReal(graph.getNodes()[input].getName(), res[i]))
      {
        delete[] res;
        return oms_status_error;
      }
    }

    // calculate residuals
    maxRes = 0.0;
    double value;
    for (int i=0; i<size; ++i)
    {
      int output = SCC[i].first;
      if (oms_status_ok != getReal(graph.getNodes()[output].getName(), value))
      {
        delete[] res;
        return oms_status_error;
      }
      res[i] -= value;

      if (fabs(res[i]) > maxRes)
        maxRes = fabs(res[i]);
    }
  } while(maxRes > tolerance && it < maxIterations);

  delete[] res;

  if (it >= maxIterations)
    return logError("max. number of iterations (" + std::to_string(maxIterations) + ") exceeded at time = " + std::to_string(getTime()));
  logDebug("CompositeModel::solveAlgLoop: maxRes: " + std::to_string(maxRes) + ", iterations: " + std::to_string(it) + " at time = " + std::to_string(getTime()));
  return oms_status_ok;
}
