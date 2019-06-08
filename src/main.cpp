// Copyright (C) 2017-2019 Rahul Kumar
// Basque center for Applied Mathematics
// This is the main file for implementing the numerical methods
// This file is part of Numerical simulation of fouth order total variation flow using IPDG methods.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
// First published : 06/04/2019

#include <dolfin.h>
#include "Tvf4.h"
#include "TVF.h" // Change
#include "inputoutput.h" // Change
#include "converge.h"  
#include "TvfGrad.h"   // Change
#include "ErrorIPDG.h"
//#include "ErrorIPDG1.h"
//#include "ErrorIPDG2.h"
#include <fstream>
#include <csignal>
#include <iostream>
#include <math.h>


#include <boost/filesystem.hpp>

using namespace dolfin;


#define Length          1.0 
#define NumberOfNodes 128
#define MaximumTimeSteps 100
#define TimeStepSize  1.0e-03

#define ALPHA 113.0
#define H 1.0/64.0
#define K 0.5


// Signal handler to throw unexpected signals 
volatile std::sig_atomic_t gSignalStatus; 
void signal_handler(int signal)
{
  gSignalStatus = signal;
}

#define TVF_SIM  // Change

class Simulation   
{
  public:
  Simulation(const std::vector<size_t> &nodeVec, struct initParams &params)
      : numMeshes(nodeVec.size())
  {
    for(size_t i=0; i<numMeshes; i++)
    {
      params.numNodes = nodeVec[i];
  #if defined(TVF_SIM) // Change
      meshes.push_back(new TVF(params));  
  #endif    
    }
    currentTime = 0.0;
    currentTimeStep = meshes[0]->currentTimeStep;
  }
  ~Simulation()
  {
    for(size_t i=0; i<numMeshes; i++)
    {
      delete meshes[i];
    }
  }
  
 void saveConvergence(std::string fbase)
  {
    std::stringstream sstream;
    std::ofstream logFile;
    std::string fpath;

 // sstream << L2norm(1-2)<<","<< L2norm(2-3)<<","<< H1norm(1-2)<<","<< H1norm(2-3) "," << DGnorm(1-2)<<","<< DGnorm(2-3);
    sstream << std::endl;

    for(int i=0; i<eL2u.size(); i++) // // Change
    {  

      sstream << eL2u[i].first<<","<< eL2u[i].second<<","<< eH1u[i].first<<","<< eH1u[i].second<<","<< eDG[i].first<<","<< eDG[i].second;
      sstream << std::endl;
    }

    fpath.assign(fbase + "convergence.txt");
      logFile.open(fpath, std::ios::trunc);
        logFile << sstream.str();
      logFile.flush();
      logFile.close();
      std::cout<<"Convergence data written to file..."<<std::endl;
  }
  void computeConvergence()
  {
    // Assigninig functions space, billinear form and linear form
    auto Vconv = std::make_shared<converge::FunctionSpace>(meshes[0]->mesh);
    auto L_h1 = std::make_shared<converge::LinearForm>(Vconv);
    auto a_h1 = std::make_shared<converge::BilinearForm>(Vconv,Vconv);
    auto Vconv2 = std::make_shared<ErrorIPDG::FunctionSpace>(meshes[0]->mesh);
    auto L_h2 = std::make_shared<ErrorIPDG::LinearForm>(Vconv2);
    auto a_h2 = std::make_shared<ErrorIPDG::BilinearForm>(Vconv2,Vconv2);
   // auto L_h3 = std::make_shared<ErrorIPDG2::LinearForm>(Vconv2);
   // auto a_h3 = std::make_shared<ErrorIPDG2::BilinearForm>(Vconv2,Vconv2);
    auto k = std::make_shared<Constant>(K);
    auto h = std::make_shared<Constant>(H);

    auto e1 = std::make_shared<dolfin::Function>(Vconv);
    auto e2 = std::make_shared<dolfin::Function>(Vconv);
    auto h1 = std::make_shared<dolfin::Function>(Vconv);
 //   auto E1 = std::make_shared<dolfin::Function>(Vconv2);
   // auto E2 = std::make_shared<dolfin::Function>(Vconv2);
    auto IPDG1 = std::make_shared<dolfin::Function>(Vconv2);
    auto IPDG2 = std::make_shared<dolfin::Function>(Vconv2);
  //  auto IPDG3 = std::make_shared<dolfin::Function>(Vconv2);
  //  auto IPDG4 = std::make_shared<dolfin::Function>(Vconv2);
    auto DG1 = std::make_shared<dolfin::Function>(Vconv2);
  //  auto DG2 = std::make_shared<dolfin::Function>(Vconv2);

    L_h1->u1      = e1;
    L_h1->u2      = e2;

    L_h2->u1      = IPDG1;
    L_h2->u2      = IPDG2;
    L_h2->k       = k;
    L_h2->h       = h;

  //  L_h3->u1      = IPDG3;
   // L_h3->u2      = IPDG4;
  //  L_h3->k       = k;
  //  L_h3->h       = h;

  // # Calculating the L2 norm
 // # Interpolate functions into (higher order) finite element space
    e1->interpolate(*meshes[0]->u);
    e2->interpolate(*meshes[1]->u);
    e1->vector()->axpy(-1.0, *(e2->vector()));
    auto result1 = e1->vector()->norm("l2");
    std::cout<<"CONVERGENCE RESULT U NORM L2 (1-2): "<<result1<<std::endl;
    e1->interpolate(*meshes[2]->u);
    e1->vector()->axpy(-1.0, *(e2->vector()));
    auto result2 = e1->vector()->norm("l2");
    std::cout<<"CONVERGENCE RESULT U NORM L2(2-3): "<<result2<<std::endl;    
    eL2u.push_back(std::make_pair(result1, result2));

// Calculating the H1 norm
    e1->interpolate(*meshes[0]->u);
    e2->interpolate(*meshes[1]->u);
    solve(*a_h1 == *L_h1, *h1);
    result1 = h1->vector()->norm("l2");
    std::cout<<"CONVERGENCE RESULT U NORM H1 (1-2): "<<result1<<std::endl;
    e1->interpolate(*meshes[2]->u);
    solve(*a_h1 == *L_h1, *h1);
    result2 = h1->vector()->norm("l2");
    std::cout<<"CONVERGENCE RESULT U NORM H1 (2-3): "<<result2<<std::endl;
    eH1u.push_back(std::make_pair(result1, result2));


// Calculating the Mesh dependent norm
 //   E1->interpolate(*meshes[0]->u);
  //  E2->interpolate(*meshes[1]->u);
//    E1->vector()->axpy(-1.0, *(E2->vector()));
    IPDG1->interpolate(*meshes[0]->u);
    IPDG2->interpolate(*meshes[1]->u);
 //   IPDG3->interpolate(*meshes[0]->u);
 //   IPDG4->interpolate(*meshes[1]->u);
    solve(*a_h2 == *L_h2, *DG1);
 //   solve(*a_h3 == *L_h3, *DG2);
    result1 =  DG1->vector()->norm("l2"); //+ DG2->vector()->norm("l2");
    std::cout<<"CONVERGENCE RESULT Mesh Dependent (1-2): "<<result1<<std::endl;
 //   E1->interpolate(*meshes[2]->u);
 //   IPDG3->interpolate(*meshes[2]->u);
    solve(*a_h2 == *L_h2, *DG1);
  //  solve(*a_h3 == *L_h3, *DG2);
    result2 = DG1->vector()->norm("l2");// + (1/pow(2,1+K))*DG2->vector()->norm("l2");
    std::cout<<"CONVERGENCE RESULT Mesh Dependent (2-3): "<<result2<<std::endl;
    eDG.push_back(std::make_pair(result1, result2));

  }
  
  std::shared_ptr<dolfin::Function> u(){return meshes[0]->u;}

  void copySolutions()
  {
    for(size_t i=0; i<numMeshes; i++)
    {
       meshes[i]->copySolutions();
    }
  }
  std::tuple<std::size_t, bool, double, double> step_u()
  {
    auto uConverged = false;
    auto uNewtonSteps = 1;

    double u_delta, tscaleU;
    for(size_t i=0; i<numMeshes; i++)
    {
       uResult = meshes[i]->step_u();
        uConverged    = get<1>(uResult);
        if(!uConverged) return uResult;
        if(get<0>(uResult) > uNewtonSteps) uNewtonSteps  = get<0>(uResult);
    }
    return std::make_tuple(uNewtonSteps, true, 0.0, 0.0);
  }
 
  void changeTimeStep(double scale)
  {
    for(size_t i=0; i<numMeshes; i++) 
    {
      meshes[i]->changeTimeStep(scale);
    }
    currentTimeStep = meshes[0]->currentTimeStep;
  }
  void resetSolutions()
  {
    for(size_t i=0; i<numMeshes; i++) 
      meshes[i]->resetSolutions();
  }
  void incCurrentTime() {currentTime += currentTimeStep;}
  double currentTime, currentTimeStep;
  private:
  #if defined(TVF_SIM)
    std::vector<TVF *> meshes;// Change 
  #endif    
  size_t numMeshes;
  std::tuple<std::size_t, bool, double, double>  uResult;
  
  std::vector<std::pair<double, double>> eL2u, eH1u, eDG;
};
////////////////////////////////////////////////////////////////////////////////
//			                       MAIN()
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  std::cout<<"ENTERED MAIN..."<<std::endl;

  // Install a signal handler
  std::signal(SIGINT, signal_handler);
  std::signal(SIGTERM, signal_handler);
           
  auto fileIO = inputoutput();

  auto test = fileIO.parseInputLine(argc, argv);
  if(0 == test) return 0;


  auto abortPath = "../abort.txt";
  auto p = boost::filesystem::path(abortPath);   // p reads clearer than argv[1] in the following code

  if (boost::filesystem::exists(p))
  {
    std::cout<<"file: "<<abortPath<<" exists. Deleting this file..."<<std::endl;
    boost::filesystem::remove(p);
  }

//**************************************************************************//
//			             CLASS OBJECTS	
//**************************************************************************//
//{
  
  //for MPI:
  std::string old_ghost_mode = parameters["ghost_mode"];
  parameters["ghost_mode"] = "shared_facet";

  std::cout << "dolfin_version(): " << dolfin_version() << std::endl;
  std::cout << "ufc_signature(): " << ufc_signature() << std::endl;

//}		
//**************************************************************************//
//			             PARAMETERS:
//**************************************************************************//
//{
  struct initParams params;

  params.gridSizeMicrons    = std::array<double,2>
                                {Length,Length};
  params.numNodes = NumberOfNodes;

  params.Alpha = ALPHA;

  params.Delta = 4.5e-2;
  params.Beta = 1.0e-6;
  params.ScaledTimeStepSize = TimeStepSize;


//}
////////////////////////////////////////////////////////////////////////////////
//                SLURM ARRAY INDEXING
////////////////////////////////////////////////////////////////////////////////
    #define NUM_ARRAY_VALS  4
    double dataVector[NUM_ARRAY_VALS] = {2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3};

    if( (fileIO.isOpuntiaCluster || fileIO.isAugsbergCluster) 
      && (fileIO.slurmArrayIndex < NUM_ARRAY_VALS) 
      && (fileIO.slurmArrayIndex >= 0) )
    {
      // do nothing
    } 
////////////////////////////////////////////////////////////////////////////////
//                FENICS VARIABLES
////////////////////////////////////////////////////////////////////////////////
//{
  std::cout<<"creating simulation object...\n";
  
  // TVF Parameters (params); Different mesh sizes
  std::vector<size_t> meshSizes = {127,31,63};
  //std::vector<size_t> meshSizes = {127,15,31};
  Simulation tvf(meshSizes, params);

  std::cout<<"created simulation object...\n";

//}
//**************************************************************************//
//   					    INITIALIZE DATA-TO-FILE:
//**************************************************************************//  
//{

  auto fbase = fileIO.initOutputFiles(params); 
  
  std::string fpath;
  fpath.assign(fbase + "tvf_sol.pvd");  dolfin::File fileu(fpath, "compressed");

  //DUMP INITIAL DATA TO FILES:
  // filet << *TVF.theta;

  fileu << *tvf.u();
//}
//**************************************************************************//
//                              MAIN LOOP
//**************************************************************************//
//{
  std::tuple<std::size_t, bool, double, double>  uResult;

  auto uConverged = false;
  auto uNewtonSteps = 1;

  double u_delta, tscaleU;

  bool changeTime         = false;

  double timingData[2*MaximumTimeSteps][5];//time, steps, currentTimeSetp, dx0, u1
  int failedConvergence = 0;


  auto TimeStep = 0;
  #define maxTimeStep 1000.0*TimeStepSize


    while (TimeStep < MaximumTimeSteps)
    {

      if (boost::filesystem::exists(p))
      {
        std::cout<<"file: "<<abortPath<<" exists...terminating..."<<std::endl;
        break;
      }

      if( (SIGINT == gSignalStatus) || (SIGTERM == gSignalStatus) )
      {
          std::cout << "SignalValue: " << gSignalStatus << '\n';
          break; //break while loop to exit sim
      }

      //new code:
    
     
      tvf.copySolutions(); // Copy the old solution 
      uResult = tvf.step_u(); // Continue to next time step
      uConverged    = get<1>(uResult);
      uNewtonSteps  = get<0>(uResult);

      if(false == uConverged)
      {//DID NOT CONVERGE:
        failedConvergence++;
        changeTime = true;
        //tscaleU = get<3>(uResult);
        std::cout<<"FAILED TO CONVERGE U (count="<<failedConvergence<<")"<<std::endl;
        std::cout<<"PREDICTED NEW TIME SCALING (THETA): "<< tscaleU<<std::endl;
        //over-ride the value to this:
        tscaleU = 0.93;
        // tscaleTheta = 0.5;
      }
      else
      {//CONVERGED:
        if(uNewtonSteps <= 2)//if converged in 1 step: scale theta 2x (potential deadlock otherwise)
        // if(0)//if converged in 1 step: scale theta 2x (potential deadlock otherwise)
        {
          changeTime = true;
          tscaleU = 1.1;
          // tscaleTheta = 1.02;
          std::cout<<"1 to 2 step CONVERGENCE, SCALING TIMESTEP: "<<tscaleU<<"x..."<<std::endl;
        }
        else
        {//converged with > 2 steps for u:
        
          changeTime = false;
          //predictor-corrector computation (now not used:)
/*          #define CORRECTION_SCALE 0.15
          tscaleU = CORRECTION_SCALE
                           * ( (sqrt(2.0)-1.0) * get<2>(uResult) )
                           / ( 2.0 * get<3>(uResult) * u_delta ); */ 
        }

      

      if(uConverged)
        {
          //TODO:  FIX THE CONV. STEPS:
          fileIO.logTimingData(tvf.currentTime, double(get<0>(uResult)),
                                tvf.currentTimeStep, get<2>(uResult), get<3>(uResult),
                                timingData);
        }
        else
        {
          changeTime = true;
          tscaleU = 0.9;
          // maxTimeStep = polyc.currentTimeStep;
        }

      }

     
      if(changeTime)
      {//IF NOT CONVERGED, RESET AND CONTINUE...
        
        tvf.changeTimeStep(tscaleU);
        if(tscaleU > 1.0)
        {
          if (tscaleU * tvf.currentTimeStep > maxTimeStep)
          {
            tscaleU = 1.0;
          }
          else
          {
            std::cout<<"CONVERGENT TIME STEP INCREASE!"<<std::endl;
          }

        }

     

        if(uConverged)

           {//CONVERGED:
             std::cout<<"CURRENT TIME (STEP): "<<tvf.currentTime<<" ("<<tvf.currentTimeStep<<")"<<std::endl;
           }
        else
           {//FAILED TO CONVERGE, RESTART TIME STEP:
          //we did not perform phi soln. if did not converge; reset solution from buffer:
             tvf.resetSolutions();
             std::cout<<"RESTARTING TIME STEP:"<<TimeStep<<std::endl;
             std::cout<<"CURRENT TIME (STEP): "<<tvf.currentTime<<" ("<<tvf.currentTimeStep<<")"<<std::endl;
             continue;//jump back to start of while loop
           }
       }
      // #RK: this is the integer counter for total iterations:
      TimeStep += 1;
      tvf.incCurrentTime(); 
      tvf.computeConvergence();

      //DECIMATION OF FILE OUTPUT:
      // if( ((0 == TimeStep%10) && fileIO.isLocalComputer )
      if( ((0 == TimeStep%10) || fileIO.isLocalComputer )
       && (0 == TimeStep%10) )//decimate  x100
      // if(0 == TimeStep%10)//decimate x10
  		{//FILE OUTPUT:
       
           fileu << *tvf.u();  
  		}

      std::cout<<"End Iteration: "<<TimeStep
      		<<" of: "<<MaximumTimeSteps
          <<" failedConvergence  = "<<failedConvergence
          <<std::endl<<"CURRENT TIME (STEP): "<<tvf.currentTime<<" ("<<tvf.currentTimeStep<<")"          
      		<<std::endl;  

    }//end while()
//}
//**************************************************************************//
//   				                   END MAIN LOOP
//**************************************************************************//
    fileIO.finalizeLogFile(timingData);
    tvf.saveConvergence(fbase);
// #endif    

}//end main()

