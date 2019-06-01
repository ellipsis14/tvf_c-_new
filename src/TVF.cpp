
#include "Tvf4.h"
#include "TVF.h"

#include "inputoutput.h"
#include "converge.h"
#include "ErrorIPDG1.h"
#include "ErrorIPDG2.h"
#include <iostream>

//NEWTON SOLVER FILES RE-WRITTEN:
/*#include "./jdolfin/jNonlinearVariationalSolver.cpp"
#include "./jdolfin/jNewtonSolver.cpp"*/

#define Length 			    1.0 

///////////////////////////////////////////////////////////////////////////////
#define NumberOfNodes	62
#define MaximumTimeSteps 1200
///////////////////////////////////////////////////////////////////////////////


// #define TimeStepSize  5.0e-03
#define TimeStepSize  1.5e-03

#define ALPHA	150.0

////////////////////////////////////////////////////////////////////////////////


/* class InitialConditions : public Expression
 {
     public:
     void eval(Array<double>& values,  const Array<double>& x) const
     {
         std::pair<double,double> range1 = {0.2,0.8};
         std::pair<double,double> range2 = {0.4,0.6};
         std::pair<double,double> range3 = {0.6,0.8};

         if(between(x[0],range1) )
         {
           if (between(x[1],range2))  
                values[0]=4.0*pow(10,-2);  //1
           else
                values[0]=-3.0*pow(10,-2);//3
         }
         else if( between(x[0],range2) )
         {
           if (between(x[1],range1))  
                values[0]=4.0*pow(10,-2);//2
           else
                values[0]=-3.0*pow(10,-2);//3
         }
      }

 }; */


 class InitialConditions : public Expression
{
public:
   void eval(Array<double>& values, const Array<double>& x) const
     {
         values[0]= x[0]*(x[0]-1.0)*(x[0]-0.5)*x[1]*(x[1]-1.0)*(x[1]-0.5)-1.0/32;
     }
};


// /////////////////////////////////////////////////////////////////////////////
//       Example-3
////////////////////////////////////////////////////////////////////////////////
// class InitialConditions : public Expression
// {
//     public:
//     void eval(Array<double>& values,  const Array<double>& x) const
//     {
//         std::pair<double,double> range1 = {0.0,0.5};
//         std::pair<double,double> range2 = {0.5,1.0};

//         if(between(x[0],range1) )
//         {
//           if (between(x[1],range1))  
//                values[0]=9*x[0]*(0.5-x[0])*(1.0-x[0])*x[1]*(1.0-x[1])*(0.5-x[1])-1.0/64;//1
//           else
//                values[0]=9*x[0]*(0.5-x[0])*(1.0-x[0])*x[1]*(1.0-x[1])*(x[1]-0.5)-1.0/64;//3
//         }
//         else if( between(x[0],range2) )
//         {
//           if (between(x[1],range1))  
//                values[0]=9*x[0]*(x[0]-0.5)*(1.0-x[0])*x[1]*(1.0-x[1])*(0.5-x[1])-1.0/64;//2
//           else
//                values[0]=9*x[0]*(1.0-x[0])*(x[0]-0.5)*x[1]*(1.0-x[1])*(x[1]-0.5)-1.0/64;//4
//         }
//      }
// };  
// /////////////////////////////////////////////////////////////////////////////
//       Example-4
////////////////////////////////////////////////////////////////////////////////
 // class InitialConditions : public Expression
 // {
 //     public:
 //     void eval(Array<double>& values,  const Array<double>& x) const
 //     {
 //         std::pair<double,double> range1 = {0.2,0.8};
 //         std::pair<double,double> range2 = {0.4,0.6};
 //         std::pair<double,double> range3 = {0.6,0.8};

 //         if(between(x[0],range1) )
 //         {
 //           if (between(x[1],range2))  
 //                values[0]=4.0*pow(10,-2);  //1
 //           else
 //                values[0]=-3.0*pow(10,-2); //3
 //         }
 //         else if( between(x[0],range2) )
 //         {
 //           if (between(x[1],range1))  
 //                values[0]=4.0*pow(10,-2);  //2
 //           else
 //                values[0]=-3.0*pow(10,-2); //3
 //         }
 //      }
 // }; 
// /////////////////////////////////////////////////////////////////////////////
//       Example-5 
///////////////////////////////////////////////////////////////////////////////////
/* class InitialConditions : public Expression
 {
     public:
     void eval(Array<double>& values,  const Array<double>& x) const
     {
         std::pair<double,double> range1 = {0.2,0.4};
         std::pair<double,double> range2 = {0.6,0.8};

         if(between(x[0],range1) )
         {
           if (between(x[1],range1))  
                values[0]=5.5*pow(10,-3);  //1
           else
                values[0]=-1.5*pow(10,-2);//5
         }
         else if( between(x[0],range2) )
         {
           if (between(x[1],range1))  
                values[0]=5.5*pow(10,-3);//2
           else
                values[0]=-1.5*pow(10,-2);//5
         }
         else if( between(x[0],range1) )
         {
           if (between(x[1],range2))  
                values[0]=5.5*pow(10,-3);//3
           else
                values[0]=-1.5*pow(10,-2);//5
         }
         else if( between(x[0],range2) )
         {
           if (between(x[1],range2))  
                values[0]=5.5*pow(10,-3);//4
           else
                values[0]=-1.5*pow(10,-2);//5
         }
      }
 };  */


 ////////////////////////////////////////////////////////////////////////////////
 // PERIODIC BOUNDARY CONDITIONS
 ///////////////////////////////////////////////////////////////////////////////
class PeriodicBC : public SubDomain
{
  public:
  // Left boundary is "target domain" G
  bool inside(const Array<double>& x, bool on_boundary) const
  { 
    return on_boundary && ( (std::abs(x[0]) < DOLFIN_EPS) || (std::abs(x[1]) < DOLFIN_EPS) ); 
  }
  // Map right boundary (H) to left boundary (G)
  void map(const Array<double>& x, Array<double>& y) const
  {
      if (x[1] > (Length - DOLFIN_EPS))
      {
          y[0] = x[0];
          y[1] = x[1] - Length;
      }
      else
      {
          y[0] = x[0] - Length;
          y[1] = x[1];
      }
  }
};

TVF::~TVF() 
{ 

}


TVF::TVF(struct initParams &params)
      : gridHeightMicrons(params.gridSizeMicrons[0]), 
        gridWidthMicrons(params.gridSizeMicrons[1]),
          numMeshNodes(params.numNodes)
{

  auto p0 = dolfin::Point(0.0,0.0);
  auto p1 = dolfin::Point(Length,Length);

  //future version for mesh creation:
  // mesh = dolfin::RectangleMesh::create({p0, p1}, {nodesW, nodesH},
  //      dolfin::CellType::Type::triangle, "right");
  mesh = std::make_shared<dolfin::RectangleMesh>(p0, p1, numMeshNodes, numMeshNodes, "crossed");

  //PERIODIC BC INITS:
  auto PBC = std::make_shared<PeriodicBC>();


// tvf FUNCTION SPACES    
  V_tvf = std::make_shared<Tvf4::FunctionSpace>(mesh, PBC);
  G_tvf = std::make_shared<TvfGrad::FunctionSpace>(mesh);

//tvf FUNCTIONS:
  u = std::make_shared<dolfin::Function>(V_tvf);  u->set_allow_extrapolation(true);	
  u0 = std::make_shared<dolfin::Function>(V_tvf);   u0->set_allow_extrapolation(true); 

//tvf CLASS OBJECTS:
  F = std::make_shared<Tvf4::ResidualForm>(V_tvf);
  J = std::make_shared<Tvf4::JacobianForm>(V_tvf,V_tvf);

// Gradient functions
  gradU = std::make_shared<Function>(G_tvf); gradU->set_allow_extrapolation(true);

  //gradient CLASS OBJECTS:
  L_gradU = std::make_shared<TvfGrad::LinearForm>(G_tvf);
  a_gradU = std::make_shared<TvfGrad::BilinearForm>(G_tvf,G_tvf);  
//**************************************************************************//
//			VARIABLES:
//**************************************************************************//




// Attaching values to the function
    
    auto Alpha = std::make_shared<Constant>(params.Alpha);
    auto Beta = std::make_shared<Constant>(params.Beta);
    auto Delta = std::make_shared<Constant>(params.Delta);
    auto ScaledTimeStepSize = std::make_shared<Constant>(params.ScaledTimeStepSize);
    auto ScaledTimeStepSize2x   = std::make_shared<Constant>(2.0*TimeStepSize);
    auto ScaledTimeStepSize10x   = std::make_shared<Constant>(10.0*TimeStepSize);
/////////////////////////////////////////////////////////////////////////////////////////// 
//            F : Total Variation Solution   
///////////////////////////////////////////////////////////////////////////////////////////
    F->Beta = Beta;  J->Beta = Beta;
    F->Alpha = Alpha; J->Alpha = Alpha;
    F->ScaledTimeStepSize = ScaledTimeStepSize; J->ScaledTimeStepSize = ScaledTimeStepSize;
    F->Delta =  Delta; J->Delta = Delta;
    
////////////////////////////////////////////////////////////////////////////////
//  PROBLEM/SOLVER DECLARATIONS:
////////////////////////////////////////////////////////////////////////////////

    F->u = u;  J->u = u; F->u0 = u0;

// GRADIENT EXPRESSION ASSIGNMENTS:
    L_gradU->solution = u;
// Necessary vector object to hold boundary condition (in our case its empty)
    std::vector<std::shared_ptr<const DirichletBC>> bcs;
    std::vector<std::shared_ptr<const DirichletBC>> nBC;



// # Create Nonlinear Problem and Solver
    Problem_Tvf			    = std::make_shared<NonlinearVariationalProblem>(F, u, bcs, J);
    Solver_Tvf  			    = std::make_shared<jNonlinearVariationalSolver>(Problem_Tvf);

    Solver_Tvf->parameters("newton_solver")["error_on_nonconvergence"]         = false; 

    Solver_Tvf->parameters("newton_solver")["maximum_iterations"]              = 150; 


    //**************************************************************************//
//                INITIAL CONDITIONS
//**************************************************************************//
  InitialConditions u_initial;
  u->interpolate(u_initial);

  currentTimeStep  = *ScaledTimeStepSize;

  guLVP = std::make_shared<dolfin::LinearVariationalProblem>( a_gradU, L_gradU, gradU, nBC);
  guLVS = std::make_shared<dolfin::LinearVariationalSolver>(guLVP);
  guLVS->solve();



  //extra buffers of prev. soln. and differences:
  u_prev    = std::make_shared<dolfin::Function>(V_tvf);    u_prev->set_allow_extrapolation(true);
  uDiff   = std::make_shared<dolfin::Function>(V_tvf);    uDiff->set_allow_extrapolation(true);

}

void 
TVF::changeTimeStep(double scale)
{
    currentTimeStep *= scale;
    newScaledTimeStepSize.reset( new Constant(currentTimeStep)) ;
    F->ScaledTimeStepSize = newScaledTimeStepSize;  J->ScaledTimeStepSize = newScaledTimeStepSize;
    std::cout<<"CHANGING TIME STEP TO: "<<currentTimeStep<<std::endl;
}

void 
TVF::resetSolutions()
{
  std::vector<double> solution_vector;
  u_prev->vector()->get_local(solution_vector);
  u->vector()->set_local(solution_vector);
  solution_vector.clear(); 
}
void 
TVF::copySolutions()
{
  std::vector<double> solution_vector;
        //save tho old solution

   u->vector()->get_local(solution_vector);
   uDiff->vector()->set_local(solution_vector);
   u_prev->vector()->set_local(solution_vector);
   u0->vector()->set_local(solution_vector); 
   solution_vector.clear(); 

}

std::tuple<std::size_t, bool, double, double>
TVF::step_u()
{        
     // continue it to next time step:
      std::tuple<std::size_t, bool, double, double>  uResult;
      std::cout<<"SolverU->solve();";  
      uResult = Solver_Tvf->solve();
      auto uConverged    = get<1>(uResult);
      if(uConverged)
      {
        std::cout<<"solve(a_U == L_U, gradU); ";  
        guLVS->solve();
      }

      return uResult;

}



std::string 
TVF::printMeshData(void)
{

  sstream.str("");

  sstream<<"mesh->num_cells(): "<<mesh->num_cells()<<std::endl;
  sstream<<"mesh->num_vertices(): "<<mesh->num_vertices()<<std::endl;
  sstream<<"mesh->num_edges(): "<<mesh->num_edges()<<std::endl;
  sstream<<"mesh->num_faces(): "<<mesh->num_faces()<<std::endl;
  sstream<<"mesh->num_facets(): "<<mesh->num_facets()<<std::endl;


  return sstream.str();
}

#if 0
//**************************************************************************
//                PREDICTOR-CORRECTOR ALGORITHM
//**************************************************************************//
//{
  std::tuple<std::size_t, bool, double, double>  thetaResult;
  std::tuple<std::size_t, bool, double, double>  concResult;
  std::tuple<std::size_t, bool, double, double>  phiResult;

  //extra buffers of prev. soln. and differences:
  auto theta_0    = std::make_shared<Function>(V_Theta);  theta_0->set_allow_extrapolation(true);
  auto c_0      = std::make_shared<Function>(V_Conc);     c_0->set_allow_extrapolation(true);
  auto phi_0      = std::make_shared<Function>(V_Phi);    phi_0->set_allow_extrapolation(true);
  auto thetaDiff  = std::make_shared<Function>(V_Theta);  thetaDiff->set_allow_extrapolation(true);
  auto cDiff    = std::make_shared<Function>(V_Conc);     cDiff->set_allow_extrapolation(true);
  auto phiDiff    = std::make_shared<Function>(V_Phi);    phiDiff->set_allow_extrapolation(true);

  double phi_delta, theta_delta, tscalePhi, tscaleTheta;

  bool changeTime         = false;
  double currentTimeStep  = *ScaledTimeStepSize;
  double currentTime      = 0.0;
  double prevStep         = 0.0;
  std::cout<<"INITIAL TIME STEP: "<<currentTimeStep<<std::endl;

  double timingData[2*MaximumTimeSteps][5];//time, steps, currentTimeSetp, dx0, theta1
  int failedConvergence = 0;


//}
#endif   


                           

