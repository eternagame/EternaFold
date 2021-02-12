//////////////////////////////////////////////////////////////////////
// BundleMethod.ipp
//
// This file contains an implementation bundle method optimization 
// algorithm
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// BundleMethod::BundleMethod()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

BundleMethod::BundleMethod
(
    const int     MAX_ITERATIONS,                     // maximum number of iterations to run subgradient method
    const double lambda
) :
  MAX_ITERATIONS(MAX_ITERATIONS),
  lambda(lambda)
{
}


//////////////////////////////////////////////////////////////////////
// BundleMethod::Minimize()
//
// Implementation of bundle methods for optimization
//////////////////////////////////////////////////////////////////////

RealT BundleMethod::Minimize(std::vector<RealT> &x)
{

  CBMRMInnerSolver* innerSolver;  // pointer to inner solver object 
#ifdef DAIFLETCHER              
  cout << "Using daifletcher as inner solver\n";
  innerSolver = new CL2N2_DaiFletcherPGM(lambda);
#else
  cout << "Using loqo as inner solver\n";
  innerSolver = new CL2N2_prLOQO(lambda);
#endif
  int verbosity        = 2;
  int convergenceLog   = 0;
  unsigned int maxNumOfIter     = 10000;
  double epsilonTol       = 1e-4;
  double gammaTol         = 0.0;
  std::string checkpointPrefix = "model.checkpoint";
  //unsigned int checkpointInterval = 1000000;  // no checkpoint by default
  //unsigned int checkpointMode = 2;

  std::vector<RealT> g;
  ComputeSubgradient(g, x);
  RealT f = ComputeFunction(x);
  

  /* Start from here */
  Scalar temp;
  /* Copy g to the matrix format */
  const int col = g.size();
  const int row = 1;
  TheMatrix a(row, col, SML::DENSE);   // gradient vector                                                                         
  TheMatrix w(row, col, SML::DENSE);   // weight vector                                                                         
  TheMatrix w_best(row, col, SML::DENSE);   // best weight vector                                                                         
  a.Zero();
  w.Zero();
  w_best.Zero();
  /*for (int i=0;i<col; i++){
    x[i] = i;
    }*/
  for (int i=0;i<col; i++){
    w.Set(i, x[i]);
  }

  /*cout <<"print X init\n";
  for (int i=0;i<col; i++){
    cout << x[i] << " ";
  }
  cout << "\n";
  cout <<"print W init\n";
  for (int i=0;i<col; i++){
    w.Get(i, temp);
    cout << temp << " ";
  }
  cout << "\n";*/

  // Timers (CPU and wall-clock)
  CTimer totalTime;             // total runtime of the training
  CTimer innerSolverTime;       // time for inner optimization (e.g., QP or LP)
  CTimer lossAndGradientTime;   // time for loss and gradient computation
  
  unsigned int iter = 0;              // iteration count
  Scalar loss = 0.0;            // loss function value        
  Scalar exactObjVal = 0.0;           // (exact) objective function value
  Scalar approxObjVal = -SML::INFTY;    // convex lower-bound (approximate) of objective function value
  Scalar minExactObjVal = SML::INFTY; // minimum of all previously evaluated (exact) objective function value
  Scalar regVal = 0.0;          // value of the regularizer term e.g., 0.5*w'*w
  Scalar epsilon = 0.0;         // (duality) gap := exactObjVal - approxObjVal
  Scalar gamma = 0.0;           // := minExactObjVal - approxObjVal
  double innerSolverTol = 1.0;  // optimization tolerance for inner solver
  
  ofstream lossFp;              // keep loss values
  ofstream exactObjValFp;       // keep exactObjVal values
  ofstream approxObjValFp;      // keep approxObjVal values
  ofstream regValFp;            // keep regVal values
  ofstream epsilonFp;           // keep epsilon values
  ofstream gammaFp;             // keep gamma values
  
  // convergence log files
  if(convergenceLog) 
    {
      lossFp.open("loss.dump");    
      exactObjValFp.open("exactobj.dump");
      approxObjValFp.open("approxobj.dump");
      regValFp.open("regval.dump");
      epsilonFp.open("epsilon.dump");
      gammaFp.open("gamma.dump");	
    }
  
  // start training
  totalTime.Start();
  while(1)
    {
      /*if (iter > 1)
	return f;*/
      iter++;
      // column generation
      lossAndGradientTime.Start();      
      for (int i=0;i<col; i++){
	w.Get(i, temp);
	x[i] = temp;
	a.Get(i, temp);
	g[i] = temp;
      }
      /*cout << "col = "<< col <<"\n";
      cout << "lambda = " << lambda << "\n";
      cout << " X == W " << "\n";
      cout <<"print X, iter = " << iter << "\n";
      for (int i=0;i<col; i++){
	cout << x[i] << " ";
      }
      cout << "\n";
      cout <<"print W, iter = "<< iter <<"\n";
      for (int i=0;i<col; i++){
	w.Get(i, temp);
	cout << temp << " ";
      }
      cout << "\n";*/

      ComputeSubgradient(g, x);
      RealT f = ComputeFunction(x);
      loss = f;
      
      for (int i=0;i<col; i++){
	temp = g[i];
	a.Set(i, temp);
      }
      /*cout << " f == loss " << "\n";
      cout << " A == G " << "\n";
      cout << "f = " << f << "\n";
      cout << "loss = "<< loss << "\n";
      cout << "print G iter = "<< iter <<" \n";
      for (int i=0;i<col; i++){
	cout << g[i] << " " ;
      }
      cout << "\n";
      cout << "print A iter = " << iter <<"\n";
      for (int i=0;i<col; i++){
	a.Get(i, temp);
	cout << temp << " " ;
      }
      cout << "\n";*/
      //cout << "Finish copying\n";
      lossAndGradientTime.Stop();
      
      // update convergence monitor
      regVal = innerSolver->ComputeRegularizerValue(w);
      exactObjVal = loss + regVal;
      //minExactObjVal = std::min(minExactObjVal,exactObjVal);
      if (minExactObjVal > exactObjVal){
	w_best.Assign(w);
	minExactObjVal = exactObjVal;
      }
      epsilon = exactObjVal - approxObjVal;
      gamma = minExactObjVal - approxObjVal;
      
      // dump convergence statistics into files
      if(convergenceLog) 
	{
	  lossFp         << loss         << endl;
	  //xiFp           << xi           << endl;
	  exactObjValFp  << exactObjVal  << endl;   
	  approxObjValFp << approxObjVal << endl;
	  regValFp       << regVal       << endl;    
	  epsilonFp      << epsilon      << endl;
	  gammaFp        << gamma        << endl;
	}
      
      // dump convergence statistics on stdout
      if(verbosity < 1) 
	{
	  printf(".");
	  if(iter%100 == 0) 
	    printf("%d",iter);
	  fflush(stdout);
	}
      else if(verbosity == 1)
	printf("#%d   eps %.6e   loss %.6e   reg %.6e\n",iter, epsilon, loss, regVal);      
      else if(verbosity > 1)
	printf("#%d   f %.6e pobj %.6e   aobj %.6e   eps %.6e   gam %.6e   loss %.6e    reg %.6e\n", iter, minExactObjVal, exactObjVal, approxObjVal, epsilon, gamma, loss, regVal);   
      
      // stopping criteria
      if((iter >= 2) && ((gamma < gammaTol) || (epsilon < epsilonTol)))
	{
	  break;
	}
      if(iter >= maxNumOfIter)
	{ 
	  printf("\nWARNING: program exceeded maximum number of iterations (%d) !\n", maxNumOfIter);
	  break;
	}
      
      // adjust inner solver optimization tolerance
      innerSolverTol = std::min(innerSolverTol, std::max((double)epsilon, (double)epsilonTol));
      innerSolverTol = std::min(innerSolverTol, std::max((double)gamma, (double)gammaTol));
      innerSolver->SetTolerance(innerSolverTol*0.5);        
			
      innerSolverTime.Start();
      //innerSolver->Solve(w, a, loss, xi, regVal, approxObjVal);
      innerSolver->Solve(w, a, loss, approxObjVal);
      innerSolverTime.Stop();
    }
  
  // legends
  if(verbosity >= 1) 
    {
      printf("\nLegends::\n");
      if(verbosity > 1)
	printf("pobj: primal objective function value\naobj: approximate objective function value\n");
      //printf("eps: epsilon (approximation error) \ngam: lower bound on eps \nloss: loss function value \nxi: approximation to loss \nreg: regularizer value\n");
      printf("eps: epsilon (approximation error) \ngam: lower bound on eps \nloss: loss function value \nreg: regularizer value\n");
    }
  
  w.Assign(w_best);
  Scalar norm1 = 0, norm2 = 0, norminf = 0;
  w.Norm1(norm1);
  w.Norm2(norm2);
  w.NormInf(norminf);
  printf("\n");
  printf("No. of iterations:  %d\n",iter);
  printf("Primal obj. val.: %.6e\n",exactObjVal);
  printf("Approx obj. val.: %.6e\n",approxObjVal);
  printf("Primal - Approx.: %.6e\n",exactObjVal-approxObjVal);
  printf("Loss:             %.6e\n",loss);
  printf("|w|_1:            %.6e\n",norm1);
  printf("|w|_2:            %.6e\n",norm2);
  printf("|w|_oo:           %.6e\n",norminf);
  
  totalTime.Stop();
  // end of training
  
  // display timing profile
  printf("\nCPU seconds in:\n");
  printf("1. loss and gradient: %8.2f\n", lossAndGradientTime.CPUTotal());
  printf("2. solver:            %8.2f\n", innerSolverTime.CPUTotal()); 
  printf("               Total: %8.2f\n", totalTime.CPUTotal());
  
  // clean up
  if(convergenceLog)
    {
      lossFp.close();
      //xiFp.close();
      exactObjValFp.close();
      approxObjValFp.close();
      regValFp.close();
      epsilonFp.close();
      gammaFp.close();
    }   
  //cout << "Start copying\n";
  f = minExactObjVal;
  for (int i=1;i<col; i++){
    w_best.Get(i, temp);
    x[i] = temp;
  }
  /*cout <<"X\n";
    for (int i=1;i<col; i++){
    cout << x[i] << " ";
    }
    cout <<"\n";*/
  //cout << "Finish copying\n";

  return f;
}

