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
    logDebug("DEBUGGING: Entering VariableStep solver");
    int fmuIndex = 0;
    std::map<ComRef, Component*> FMUcomponents;
    std::map<ComRef, Component*> canGetAndSetStateFMUcomponents;
    std::map<ComRef, Component*> noneFMUcomponents;
    std::vector<double> nominalOutput;
    std::vector<double> nominalInput;
    bool firstTime = true;
    unsigned int rollbackCounter = 0;
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
        nominalInput.clear();
        nominalOutput.clear();
        firstTime = false;
        for (const auto& component : getComponents()) // Collect all FMUs in maps
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
        if(canGetAndSetStateFMUcomponents.size() == 0) {logError("If no FMUs can get/set states, Variable Step solver can't be used."); return oms_status_error;}

        const std::vector< std::vector< std::pair<int, int> > >& sortedConnections = outputsGraph.getSortedConnections();
        for(int i=0; i<sortedConnections.size(); i++) // Collect all nominal values for use later in simulation.
        {
          if (sortedConnections[i].size() == 1)
          {
            logDebug("DEBUGGING: Size of sortedConnections[i] is: "+std::to_string(sortedConnections[i].size()));
            int input = sortedConnections[i][0].second;
            oms3::ComRef inputName(outputsGraph.getNodes()[input].getName());
            oms3::ComRef inputModel = inputName.pop_front();
            logDebug(inputModel);
            int output = sortedConnections[i][0].first;
            oms3::ComRef outputName(outputsGraph.getNodes()[output].getName());
            oms3::ComRef outputModel = outputName.pop_front();
            logDebug(outputModel);
            std::map<ComRef, Component*>::iterator inIterator = canGetAndSetStateFMUcomponents.find(inputModel);
            std::map<ComRef, Component*>::iterator outIterator = canGetAndSetStateFMUcomponents.find(outputModel);
            if (inIterator != canGetAndSetStateFMUcomponents.end())
            {
              if(outIterator != canGetAndSetStateFMUcomponents.end())
              {
                logDebug("DEBUGGING: We have found a connection with 2 FMUs thats good, lets get nominals.");
                if (outputsGraph.getNodes()[input].getType() == oms_signal_type_real)
                {
                  logDebug("DEBUGGING: Input is real, lets go.");
                  double inValue = 0.0;
                  if (oms_status_ok != getReal(outputsGraph.getNodes()[input].getName(), inValue)) return oms_status_error;          
                  logDebug("DEBUGGING: We just got real, do we need to?.");          
                  fmi2_import_variable_list_t *varList = fmi2_import_get_variable_list(dynamic_cast<ComponentFMUCS*>(inIterator->second)->getFMU(), 0);
                  logDebug("DEBUGGING: We just got varList");          
                  size_t varListSize = fmi2_import_get_variable_list_size(varList);
                  logDebug("DEBUGGING: Starting for loop");          
                  for (int n = 0; n < varListSize; n++)
                  {
                    fmi2_import_variable_t* var = fmi2_import_get_variable(varList, n);
                    logDebug("DEBUGGING: Checking if it corresponds.");
                    logDebug(outputsGraph.getNodes()[input].getName());
                    logDebug("DEBUGGING: var name is: "+std::string(fmi2_import_get_variable_name(var)));
                    if (inputName == std::string(fmi2_import_get_variable_name(var)) || outputName == std::string(fmi2_import_get_variable_name(var)))
                    {
                      logDebug("DEBUGGING: Found a real nominal that corresponds.");     
                      fmi2_import_real_variable_t* varReal = fmi2_import_get_variable_as_real(var);
                      double nominal = fmi2_import_get_real_variable_nominal(varReal);
                      nominalInput.push_back(nominal);                        
                    }                      
                  }
                  logDebug("DEBUGGING: We have looped through inputs. Lets clear varlist. ");
                  fmi2_import_free_variable_list(varList);               
                  logDebug("DEBUGGING: Varlist clear.");                
                }                
                
                if (outputsGraph.getNodes()[output].getType() == oms_signal_type_real)
                {
                  logDebug("DEBUGGING: Output is real, lets go.");
                  double outValue = 0.0;
                  if (oms_status_ok != getReal(outputsGraph.getNodes()[output].getName(), outValue)) return oms_status_error;   
                  logDebug("DEBUGGING: We just got real, do we need to?.");                           
                  fmi2_import_variable_list_t *varList = fmi2_import_get_variable_list(dynamic_cast<ComponentFMUCS*>(outIterator->second)->getFMU(), 0);
                  size_t varListSize = fmi2_import_get_variable_list_size(varList);
                  logDebug("DEBUGGING: Starting for loop");     
                  for (int n = 0; n < varListSize; n++)
                  {
                    fmi2_import_variable_t* var = fmi2_import_get_variable(varList, n);
                    if (inputName == std::string(fmi2_import_get_variable_name(var)) || outputName == std::string(fmi2_import_get_variable_name(var)))
                    {
                      fmi2_import_real_variable_t* varReal = fmi2_import_get_variable_as_real(var);
                      double nominal = fmi2_import_get_real_variable_nominal(varReal);
                      nominalOutput.push_back(nominal);                                           
                    }                      
                  }
                  logDebug("DEBUGGING: We have looped through outputs. Lets clear varlist. ");
                  fmi2_import_free_variable_list(varList);               
                  logDebug("DEBUGGING: Varlist clear.");
                }
                
              }
            }
          }
          else
            return oms_status_error; //Algebraic loop. TODO fix            
        }
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

      // get inputs and outputs at the end of all steps.
      std::vector<double> inputVect;
      std::vector<double> outputVect;
      if(oms_status_ok != getInputAndOutput(outputsGraph,inputVect,outputVect,canGetAndSetStateFMUcomponents)) return oms_status_error;

      // get outputs after communication steps
      //if(oms_status_ok != getOutput(outputsGraph,outputVect,FMUcomponents)) return oms_status_error;
      /*
      if(oms_status_ok != getInput(outputsGraph,inputVect,FMUcomponents)) return oms_status_error;
      logDebug("DEBUGGING: We just filled the input vector.");
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
      if(oms_status_ok != getOutput(outputsGraph,outputVect,FMUcomponents)) return oms_status_error;*/
      
      logDebug("DEBUGGING: Lets do Error control");
      if (inputVect.size() != outputVect.size() || inputVect.size() != nominalInput.size() || nominalOutput.size() != nominalInput.size()) return oms_status_error;

      double biggestDifferance = 0.0;
      double changeTolerance = 1e-2;
      double safety_factor = 0.90;
      for (int n=0; n < inputVect.size();n++) // Calculate error in the FMUs we do error_control on.
      {
        double error;
        double nominal;
        bool hack_error = false;        
        if (nominalInput[n] == nominalOutput[n])
        {
          nominal = nominalInput[n];
        }
        else if (nominalInput[n] == 1.0)
        {
          nominal = nominalOutput[n];
        }
        else if (nominalOutput[n] == 1.0)
        {
          nominal = nominalInput[n];
        }
        else
        {
          logError("Nominal value for input and output are specified differently. Will use Output nominal, be advised this might provide unwanted behavior.");
          nominal = nominalOutput[n];
        }          
        error = fabs(inputVect[n]-outputVect[n]);
        logDebug("DEBUGGING: Error is:"+std::to_string(error)+" and Nominal is: "+std::to_string(nominal));
        if (nominal == 1.0 && hack_error)
        {
          if (error*2/(fabs(outputVect[n])+fabs(inputVect[n])) > biggestDifferance)
          {
            if (fabs(outputVect[n]+inputVect[n]) > changeTolerance)
            {
              biggestDifferance = error*2/(fabs(outputVect[n])+fabs(inputVect[n])); // scale with mean?
              logDebug("DEBUGGING: scaled error is: " + std::to_string(error*2/(fabs(outputVect[n])+fabs(inputVect[n]))) + " New biggest Differance is: " + std::to_string(biggestDifferance));
            }
            else if(error > biggestDifferance)
            {
              biggestDifferance = error;
              logDebug("DEBUGGING: error is: " + std::to_string(error) + " New Biggest Differance is: " + std::to_string(biggestDifferance));
            }
          }
        }
        else
        {
          if (error/nominal > biggestDifferance)
            biggestDifferance = error/nominal;          
          logDebug("DEBUGGING: Scaling error to nominal and editing biggestDifferance");
        }
          
        
      }
      double fixRatio = safety_factor*changeTolerance/biggestDifferance;
      logDebug("DEBUGGING: fixRatio is: " + std::to_string(fixRatio));
      if (biggestDifferance > changeTolerance) //Going to rollback.
      {
        logDebug("DEBUGGING: Rollbacking New h is: " + std::to_string(stepSize));
        //Fix fmus
        for (int i=0; i<fmiImportVect.size(); ++i) // Reset all FMU states
        {
          fmi_status = fmi2_import_set_fmu_state(fmiImportVect[i], sVect[i]);
        }
        //Fix time
        time = tNext-stepSize;
        for (const auto& component : getComponents())
        {
          if (oms_component_fmu == component.second->getType())
          {
            dynamic_cast<ComponentFMUCS*>(component.second)->setFmuTime(time);
          }
        }
        //Fix Steptime.
        stepSize = stepSize*pow(fixRatio,0.25);
        rollbackCounter++;
      }
      else // Not going to rollback.
      {
        if (safety_factor*changeTolerance/biggestDifferance > 1)
          stepSize = stepSize*pow(fixRatio,0.25);

        logDebug("DEBUGGING: Not Rollbacking, h is: " + std::to_string(stepSize));
        if (!FMUcomponents.empty())
        {
          for (const auto& component : FMUcomponents) // These FMUs cant rollback, so only simulating them when we have decided on a step to take.
          {
            status = component.second->stepUntil(tNext);
            if (oms_status_ok != status)
            {
              if (cb)
                cb(modelName.c_str(), tNext, status);
              return status;
            }
          }
        }
        time = tNext;
        if (isTopLevelSystem())
        {
          getModel()->setStepAndRollIterator(stepSize,rollbackCounter);
          getModel()->emit(time);
        }
        updateInputs(outputsGraph); // updateCanGetFMUs(outputsGraph,canGetAndSetStateFMUcomponents);
        if (isTopLevelSystem())
        {
          //logDebug("DEBUGGING: Emitting h: "+stepSize+" and rollbackCounter:"+std::to_string(rollbackCounter));
          getModel()->setStepAndRollIterator(stepSize,rollbackCounter);
          getModel()->emit(time);
        }
        rollbackCounter = 0;
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
    logDebug("DEBUGGING: Entering FixedStep solver");
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



oms_status_enu_t oms3::SystemWC::getInputAndOutput(oms3::DirectedGraph& graph, std::vector<double>& inputVect,std::vector<double>& outputVect,std::map<ComRef, Component*> FMUcomponents)
{
  const std::vector< std::vector< std::pair<int, int> > >& sortedConnections = graph.getSortedConnections();
  inputVect.clear();
  int inCount = 0;
  outputVect.clear();
  int outCount = 0;
    for(int i=0; i<sortedConnections.size(); i++)
    {
      if (sortedConnections[i].size() == 1)
      {
        logDebug("DEBUGGING: Size of sortedConnections[i] is: "+std::to_string(sortedConnections[i].size()));
        int input = sortedConnections[i][0].second;
        oms3::ComRef inputName(graph.getNodes()[input].getName());
        oms3::ComRef inputModel = inputName.pop_front();
        logDebug(inputModel);
        int output = sortedConnections[i][0].first;
        oms3::ComRef outputName(graph.getNodes()[output].getName());
        oms3::ComRef outputModel = outputName.pop_front();
        logDebug(outputModel);
        if (FMUcomponents.find(inputModel) != FMUcomponents.end())
        {
          if(FMUcomponents.find(outputModel) != FMUcomponents.end())
          {
            if (graph.getNodes()[input].getType() == oms_signal_type_real)
            {
              logDebug("DEBUGGING: found a real input");
              double inValue = 0.0;
              if (oms_status_ok != getReal(graph.getNodes()[input].getName(), inValue)) return oms_status_error;
              inputVect.push_back(inValue);
              inCount++;
            }
            if (graph.getNodes()[output].getType() == oms_signal_type_real)
            {
              logDebug("DEBUGGING: found a real output");
              double outValue = 0.0;
              if (oms_status_ok != getReal(graph.getNodes()[output].getName(), outValue)) return oms_status_error;
              outputVect.push_back(outValue);    
              outCount++;
            }
          }
        }
      }
      else
      {
        logDebug("DEBUGGING: Exiting cuz algebraic loop!");
        return oms_status_error; // algebraic loop: TODO
      }
      logDebug("DEBUGGING: we have added "+std::to_string(inCount)+" inputs and "+std::to_string(outCount)+" outputs to the vectors.");
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
      if (sortedConnections[i].size() == 1)
      {
        int output = sortedConnections[i][0].first;
        int input = sortedConnections[i][0].second;
        oms3::ComRef outputName(graph.getNodes()[output].getName());
        oms3::ComRef inputName(graph.getNodes()[input].getName());
        oms3::ComRef outputModel = outputName.pop_front();
        oms3::ComRef inputModel = inputName.pop_front();
        if (FMUcomponents.find(inputModel) == FMUcomponents.end())
        {
          if (FMUcomponents.find(outputModel) == FMUcomponents.end())
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
        return oms_status_error; // Algebraic loop, TODO.
    }
  }
  else
  {
    for(int i=0; i<sortedConnections.size(); i++)
    {
      if (sortedConnections[i].size() == 1)
      {
        int output = sortedConnections[i][0].first;
        if (graph.getNodes()[output].getType() == oms_signal_type_real)
        {
          double outValue = 0.0;
          if (oms_status_ok != getReal(graph.getNodes()[output].getName(), outValue)) return oms_status_error;
          outputVect.push_back(outValue);
        }
      }
      else
        return oms_status_error; // Algebraic loop, TODO.
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
      if (FMUcomponents.find(inputModel) != FMUcomponents.end() && FMUcomponents.find(outputModel) != FMUcomponents.end())
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
        if (FMUcomponents.find(inputModel) != FMUcomponents.end() || FMUcomponents.find(outputModel) != FMUcomponents.end())
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
      logDebug("DEBUGGING: Lets check which connection it is.");
      if (FMUcomponents.find(inputModel) != FMUcomponents.end() || FMUcomponents.find(outputModel) != FMUcomponents.end())
      {
        if (graph.getNodes()[input].getType() == oms_signal_type_real)
        {
          double value = 0.0;
          logDebug("DEBUGGING: We found a link between two fmus where atleast one cant handle get and set, so lets update this (real)");
          if (oms_status_ok != getReal(graph.getNodes()[output].getName(), value)) return oms_status_error;
          if (oms_status_ok != setReal(graph.getNodes()[input].getName(), value)) return oms_status_error;

          // derivatives
          if (derBuffer)
          {
            unsigned int order;
            if (oms_status_ok == getRealOutputDerivative(graph.getNodes()[output].getName(), derBuffer, order))
            {
              logDebug(graph.getNodes()[output].getName() + " -> " + graph.getNodes()[input].getName() + ": " + std::to_string(derBuffer[0]));
              if (oms_status_ok != setRealInputDerivative(graph.getNodes()[input].getName(), derBuffer, order)) return oms_status_error;
            }
          }
          logDebug("DEBUGGING: Sorted the Derivative");
        }
        else if (graph.getNodes()[input].getType() == oms_signal_type_integer)
        {
          logDebug("DEBUGGING: We found a link between two fmus where atleast one cant handle get and set, so lets update this (int)");
          int value = 0.0;
          if (oms_status_ok != getInteger(graph.getNodes()[output].getName(), value)) return oms_status_error;
          if (oms_status_ok != setInteger(graph.getNodes()[input].getName(), value)) return oms_status_error;
        }
        else if (graph.getNodes()[input].getType() == oms_signal_type_boolean)
        {
          logDebug("DEBUGGING: We found a link between two fmus where atleast one cant handle get and set, so lets update this (boolean)");
          bool value = 0.0;
          if (oms_status_ok != getBoolean(graph.getNodes()[output].getName(), value)) return oms_status_error;
          if (oms_status_ok != setBoolean(graph.getNodes()[input].getName(), value)) return oms_status_error;
        }
        else
        {
          logDebug("DEBUGGING: We found a link between two fmus where atleast one cant handle get and set, BUT THE link isnt real/int/boolean, casting error.");
          return logError_InternalError;
        }

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
