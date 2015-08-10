GLOBALS_SECTION
  #include <fvar.hpp>
  //---  Estimating individual animal movement from observation networks
  //---  Martin W. Pedersen and Kevin C. Weng
  //---  Supplementary material S2, ADMB code
  //---  08.05.2013
  //---  Script for estimating OU process with spatial HMM and observation network data
  //---  Requires ADMB to compile, which can be downloaded freely from: http://www.admb-project.org/downloads
  //---  Requires about 1.5 GB free memory (RAM) to run with provided real data set
  //---  Takes about 5 min to analyse real data on a Intel i7 CPU @ 2.7 GHz

  // Receiver detection function, sigmoidal
  double detectfunSigmoid(double A, double pmax, double D50, double D95){
    return pmax/(1.0 + exp(log(19.0)*(A-D50)/(D95-D50)));
  }
  // Inverse receiver detection function, inverse sigmoid
  double detectfunSigmoidInv(double R, double pmax, double D50, double D95){
    return D50 + (D95-D50)/log(19.0)*log(pmax/R-1.0);
  }
  // Receiver detection function, for COVARIATE data
  double dfcov(double A, double pmax, double D50, double D95, double fac, double dref){
    double Rr = detectfunSigmoid(dref,pmax,D50,D95); // optimal detection probability at reference tag distance dref
    double xE = detectfunSigmoidInv(fac*Rr,pmax,D50,D95); 
    return detectfunSigmoid(A*xE/dref,pmax,D50,D95);
  }

  // Construct 1D Transition probability matrix for marginal distribution
  dvar_matrix makekern1Dmarg(dvariable mu, dvariable G11, dvariable noisesd, double dt, dvector xinter, dvector xg, int nx){
    dvar_matrix Pxx(1,nx,1,nx);
    Pxx.initialize();
    double x0;
    dvariable zl;
    dvariable zr;
    dvariable cdfzl;
    dvariable cdfzr;
    dvariable meanx;
    for(int i=1; i<=nx; ++i){
      x0 = xg(i);                         // Midpoint between interfaces (midpoint of grid cell)
      meanx = mu + G11*(x0-mu);           // Mean in the OU process
      zl = (xinter(1)-meanx)/noisesd;     // Standard Z-value  at left interface
      cdfzl = cumd_norm(zl);              // CDF value at left interface
      for(int j=1; j<=nx; ++j){
        zr = (xinter(j+1)-meanx)/noisesd; // Standard Z-value at right interface
	cdfzr = cumd_norm(zr);            // CDF value at right interface
        Pxx(i,j) = cdfzr - cdfzl;         // Transition probability
	zl = zr;                          // After step, right is left
	cdfzl = cdfzr;                    // After step, right is left
      }
    }
    return Pxx;
  }

  // Construct 1D Transition probability matrix for conditional distribution
  dvar_matrix makekern1Dcond(dvariable mu1, dvariable mu2, dvariable G11, dvariable G12, dvariable G22, dvariable x2, dvariable rho, dvariable s1, dvariable s2, double dt, dvector xinter, dvector xg, int nx){
    dvar_matrix Pxx(1,nx,1,nx);
    Pxx.initialize();
    double x1;
    dvariable zl;
    dvariable zr;
    dvariable cdfzl;
    dvariable cdfzr;
    dvariable mean1;
    dvariable mean2;
    dvariable meanx;
    dvariable noisesd = sqrt((1-square(rho)))*s1;; // Conditional standard deviation in the x1-direction
    for(int i=1; i<=nx; ++i){
      x1 = xg(i);                                  // Midpoint between interfaces (midpoint of grid cell)
      // Marginal means of OU process in each direction
      mean1 = mu1 + G11*(x1-mu1) + G12*(x2-mu2);   // x1-direction
      mean2 = mu2 + G22*(x2-mu2) + G12*(x1-mu1);   // x2-direction
      // Conditional mean in the x1 direction
      meanx = mean1 + s1/s2*rho*(x2-mean2);
      zl = (xinter(1)-meanx)/noisesd;              // Standard Z-value  at left interface
      cdfzl = cumd_norm(zl);                       // CDF value at left interface
      for(int j=1; j<=nx; ++j){
        zr = (xinter(j+1)-meanx)/noisesd;          // Standard Z-value at right interface
	cdfzr = cumd_norm(zr);                     // CDF value at right interface
        Pxx(i,j) = cdfzr - cdfzl;                  // Transition probability
	zl = zr;                                   // After step right is left
	cdfzl = cdfzr;                             // After step right is left
      }
    }
    return Pxx;
  }

  // --- Adjoint code for matrix multiplications starting here ---
  // Multiply vector with matrix
  void dfunvec(void);
  dvar_vector funvec(const dvar_vector& grid, const dvar_matrix P){
    int c = grid.indexmin(); 
    int C = grid.indexmax(); 
    dvar_vector ret(c,C);
    dvector vret(c,C);
    dmatrix vP = value(P);
    dvector vgrid = value(grid);
    int i,k;
    for(i=c; i<=C; ++i){
      vret(i) = 0;
      for(k=c; k<=C; ++k){
        vret(i) += vgrid(k) * vP(k,i);
      }
    }
    ret = nograd_assign(vret);
    save_identifier_string("test1");
    vgrid.save_dvector_value();
    vgrid.save_dvector_position();
    vP.save_dmatrix_value();
    vP.save_dmatrix_position();
    save_int_value(c);
    save_int_value(C);
    grid.save_dvar_vector_value();
    grid.save_dvar_vector_position();
    P.save_dvar_matrix_value();
    P.save_dvar_matrix_position();
    ret.save_dvar_vector_position();
    save_identifier_string("test2");
    gradient_structure::GRAD_STACK1->set_gradient_stack(dfunvec);
    return ret;
  }

  // Multiply vector with matrix, derivative
  void dfunvec(void){
    verify_identifier_string("test2");
    dvar_vector_position ret_pos = restore_dvar_vector_position();
    dvector dret = restore_dvar_vector_derivatives(ret_pos);
    dvar_matrix_position P_pos = restore_dvar_matrix_position();
    dmatrix P = restore_dvar_matrix_value(P_pos);
    dvar_vector_position grid_pos = restore_dvar_vector_position();
    dvector grid = restore_dvar_vector_value(grid_pos);
    int C = restore_int_value();
    int c = restore_int_value();
    dmatrix_position vP_pos = restore_dmatrix_position();
    dmatrix vP = restore_dmatrix_value(vP_pos);
    dvector_position vgrid_pos = restore_dvector_position();
    dvector vgrid=restore_dvector_value(vgrid_pos);
    verify_identifier_string("test1");
    dvector dgrid(c,C);
    dgrid.initialize();
    dmatrix dP(c,C,c,C);
    dP.initialize();
    int i,k;
    for(i=C; i>=c; --i){
      for(k=C; k>=c; --k){
        dgrid(k) += dret(i) * vP(k,i);
        dP(k,i) += dret(i) * vgrid(k); 
      }
      dret(i) = 0;
    }
    dgrid.save_dvector_derivatives(grid_pos);
    dP.save_dmatrix_derivatives(P_pos);
  }

  // Multiply matrix with matrix
  void dfunmat(void);
  dvar_matrix funmat(const dvar_matrix& grid, const dvar_matrix P){
    int r = grid.rowmin(); 
    int R = grid.rowmax(); 
    int c = grid.colmin(); 
    int C = grid.colmax();
    dvar_matrix ret(r,R,c,C);
    dmatrix vret(r,R,c,C);
    dmatrix vP=value(P);
    dmatrix vgrid=value(grid);
    int j,i,k;
    for(j=r; j<=R; ++j){
      for(i=c; i<=C; ++i){
        vret(j,i) = 0;
        for(k=c; k<=C; ++k){
          vret(j,i) += vgrid(j,k) * vP(k,i);
        }
      }
    } 
    ret=nograd_assign(vret);
    save_identifier_string("test1");
    vgrid.save_dmatrix_value();
    vgrid.save_dmatrix_position();
    vP.save_dmatrix_value();
    vP.save_dmatrix_position();
    save_int_value(r);
    save_int_value(R);
    save_int_value(c);
    save_int_value(C);
    grid.save_dvar_matrix_value();
    grid.save_dvar_matrix_position();
    P.save_dvar_matrix_value();
    P.save_dvar_matrix_position();
    ret.save_dvar_matrix_position();
    save_identifier_string("test2");
    gradient_structure::GRAD_STACK1->set_gradient_stack(dfunmat);
    return ret;
  }

  // Multiply matrix with matrix, derivative
  void dfunmat(void){
    verify_identifier_string("test2");
    dvar_matrix_position ret_pos = restore_dvar_matrix_position();
    dmatrix dret=restore_dvar_matrix_derivatives(ret_pos);
    dvar_matrix_position P_pos = restore_dvar_matrix_position();
    dmatrix P = restore_dvar_matrix_value(P_pos);
    dvar_matrix_position grid_pos = restore_dvar_matrix_position();
    dmatrix grid = restore_dvar_matrix_value(grid_pos);
    int C = restore_int_value();
    int c = restore_int_value();
    int R = restore_int_value();
    int r = restore_int_value();
    dmatrix_position vP_pos = restore_dmatrix_position();
    dmatrix vP=restore_dmatrix_value(vP_pos);
    dmatrix_position vgrid_pos = restore_dmatrix_position();
    dmatrix vgrid=restore_dmatrix_value(vgrid_pos);
    verify_identifier_string("test1");
    dmatrix dgrid(r,R,c,C);
    dgrid.initialize();
    dmatrix dP(c,C,c,C);
    dP.initialize();
    int j,i,k;
    for(j=R; j>=r; --j){
      for(i=C; i>=c; --i){
        for(k=C; k>=c; --k){
          dgrid(j,k)+=dret(j,i)*vP(k,i);
          dP(k,i)+=dret(j,i)*vgrid(j,k); 
        }
        dret(j,i)=0;
      }
    }    
    dgrid.save_dmatrix_derivatives(grid_pos);
    dP.save_dmatrix_derivatives(P_pos);
  }

 
DATA_SECTION
  // SWITCH FILE
  !! ad_comm::change_datafile_name("grid.cfg");
  // Grid variables
  init_number xmin; // Minimum x-coordinate
  init_number ymin; // Minimum y-coordinate
  init_number xmax; // Maximum x-coordinate
  init_number ymax; // Maximum y-coordinate
  init_int ngx;     // Number of grid cells in x-direction
  init_int ngy;     // Number of grid cells in y-direction

  // SWITCH FILE
  !! ad_comm::change_datafile_name("data.dat");
  // Detection function parameters
  init_number pmax;     // Parameters in sigmoid function
  init_number D50;  
  init_number D95;
  init_number dref;
  // Array information
  init_number nr;       // Number of receivers
  init_vector rx(1,nr); // x-coordinates of receivers
  init_vector ry(1,nr); // y-coordinates of receivers
  // Detection data
  init_number dt;               // Time between pings
  init_number T;                // Total time at liberty
  int nt;                       // Total number of time steps
  !! nt = T/dt + 1;
  init_number nobs;             // Number of observations
  init_vector rec(1,nobs);      // ID of receivers that received pings
  init_vector rectimes(1,nobs); // Times of received pings
  init_vector cov(1,nt);        // Covariate: detection efficiency

  // Spatial grid
  vector gx(1,ngx);
  vector gy(1,ngy);
  vector xinter(1,ngx+1);
  vector yinter(1,ngy+1);
  number dx;
  !! dx = (xmax-xmin)/(ngx-1);
  number dy;
  !! dy = (ymax-ymin)/(ngy-1);
  !! for(int i=1; i<=ngx; ++i){
  !!   gx(i) = xmin + (i-1)*dx;
  !!   xinter(i) = gx(i) - dx/2.0;
  !! }
  !! xinter(ngx+1) = gx(ngx) + dx/2.0;
  !! for(int i=1; i<=ngy; ++i){
  !!   gy(i) = ymin + (i-1)*dy;
  !!   yinter(i) = gy(i) - dy/2.0;
  !! }
  !! yinter(ngy+1) = gy(ngy) + dy/2.0;

  // Calculate data likelihood
  vector timepoints(1,nt);
  !! cout << "Calculating data likelihood array..." << endl; 
  3darray lik(1,nt,1,ngy,1,ngx); // Data likelihood array
  !! lik.initialize();
  !! lik =+ 1.0;
  !! dvector recvec(1,nr);
  !! double demap;
  !! double dist;
  !! for(int t=2; t<=nt; ++t){ // Loop over time
  !!   timepoints(t) = (t-1)*dt;
  !!   recvec.initialize();
  !!   for(int tt=1; tt<=nobs; ++tt){
  !!     if(rectimes(tt)==timepoints(t)){
  !!	   recvec(rec(tt)) = 1;
  !!     }
  !!   }
  !!   for(int k=1; k<=nr; ++k){
  !!     if(recvec(k)==1){ // If a ping was heard by receiver k
  !!   	   for(int i=1; i<=ngy; ++i){
  !!         for(int j=1; j<=ngx; ++j){
  !!           dist = sqrt(square(gx(j)-rx(k)) + square(gy(i)-ry(k)));
  !!	       demap = dfcov(dist,pmax,D50,D95,cov(t),dref);
  !!           lik(t,i,j) *= cov(t)*demap;
  !! if(lik(t,i,j)<0){cout << " i: " << i << " j: " << j << " t: " << t << " cov: " << cov(t) << " demap: " << demap <<endl;}
  !!         }
  !!       }
  !!     }
  !!     if(recvec(k)==0){
  !!   	   for(int i=1; i<=ngy; ++i){
  !!         for(int j=1; j<=ngx; ++j){
  !!           dist = sqrt(square(gx(j)-rx(k)) + square(gy(i)-ry(k)));
  !!	       demap = dfcov(dist,pmax,D50,D95,cov(t),dref);
  !!           lik(t,i,j) *= (1.0 - cov(t)*demap);
  !! if(lik(t,i,j)<0){cout << " i: " << i << " j: " << j << " t: " << t << " cov: " << cov(t) << " demap: " << demap <<endl;}
  !!         }
  !!       }
  !!     }
  !!   }
  !! }
  !! cout << "Done. Proceeding to estimation..." << endl; 

  // Grids
  3darray pred(1,nt,1,ngy,1,ngx) // Predicted distribution
  !!  pred.initialize();
  3darray phi(1,nt,1,ngy,1,ngx)  // Reconstructed distribution
  !!  phi.initialize();
  3darray smoo(1,nt,1,ngy,1,ngx) // Smoothed distribution
  !!  smoo.initialize();

PARAMETER_SECTION
  // Parameters to be estimated
  init_bounded_number mux(xmin,xmax,1);     // Home range center, x coordinate
  init_bounded_number muy(ymin,ymax,1);     // Home range center, y coordinate
  init_bounded_number logbx(-14.0,-2.0,1);  // Degree of attaction to home range center
  init_bounded_number logby(-14.0,-2.0,1);  // Degree of attaction to home range center
  init_bounded_number logitcor(-3.0,3.0,2); // Degree of attaction to home range center
  init_bounded_number logsigma(-2.0,7.0,1); // Movement parameter
  matrix Px(1,ngx,1,ngx);
  3darray Py(1,ngx,1,ngy,1,ngy);
  matrix grid(1,ngy,1,ngx);
  matrix gridT(1,ngx,1,ngy);
  vector psi(2,nt);
  // Objective function
  objective_function_value jnll;

PRELIMINARY_CALCS_SECTION 
  // Initial parameter values
  mux = 0.5*(xmin+xmax);
  muy = 0.5*(ymin+ymax);
  logbx = -9.0;
  logby = -9.0;
  logitcor = 0.0;
  logsigma = 0.4;

PROCEDURE_SECTION
  // Initialise joint negative log likelihood
  jnll = 0.0;
  // Transform parameters
  dvariable sigma = exp(logsigma);
  dvariable bx = exp(logbx);
  dvariable by = exp(logby);
  dvariable cor = (1.0*exp(logitcor)-1.0) / (1+exp(logitcor)); // Generalised logit bounds: low=-1.0, up=1.0
  dvariable bxy = -cor*sqrt(bx*by); // cor is the correlation parameter of the stationary distribution
  dvar_matrix B(1,2,1,2);
  B(1,1) = bx;
  B(2,2) = by;
  B(1,2) = bxy;
  B(2,1) = bxy;
  dvar_matrix Gamma = expm(-B*dt);
  dmatrix I = identity_matrix(1,2);
  dvar_matrix S = 0.5*square(sigma)*inv(B)*(I - Gamma*Gamma); // Covariance of conditional distribution of OU process
  dvariable rho = S(1,2)/sqrt(S(1,1)*S(2,2)); // Correlation parameter in conditional distribution of OU process

  Px = makekern1Dmarg(mux,Gamma(1,1),sqrt(S(1,1)),dt,xinter,gx,ngx); // Construct transition probability matrix for x-direction (marginal for x)
  // Construct transition probability matrices for y-direction (for y conditional on x)
  for(int j=1; j<=ngx; ++j){
    Py(j) = makekern1Dcond(muy,mux,Gamma(2,2),Gamma(1,2),Gamma(1,1),gx(j),rho,sqrt(S(2,2)),sqrt(S(1,1)),dt,yinter,gy,ngy);
  }

  // Initialize filter
  grid = lik(1) / sum(lik(1));
  pred(1) = value(grid);
  phi(1) = value(grid);
  
  dvariable lpart;
  // Run filter
  for(int t=2; t<=nt; ++t){
    gridT = trans(grid);
    for(int j=1; j<=ngx; ++j){
      gridT(j) = funvec(gridT(j),Py(j)); // Forward probabilities in y direction
    }
    grid = trans(gridT);
    grid = funmat(grid,Px);        // Forward probabilities in x direction
    pred(t) = value(grid);         // Store prediction (only value not derivative)
    grid = elem_prod(grid,lik(t)); // Data update
    lpart = sum(grid);             // Calculate likelihood contribution
    grid /= lpart+1.0e-15;         // Normalise distribution
    phi(t) = value(grid);          // Save distribution
    psi(t) = -log(lpart);          // Store likelihood contribution
  }
  jnll = sum(psi); // Calculate joint negative loglikelihood

  cout << "Negative log-likelihood value: " << jnll << " --- Pars: mux: " << mux << " muy: " << muy << " logbx: " << logbx << " logby: " << logby << " logitcor: " << logitcor << " logsigma: " << logsigma << endl;

REPORT_SECTION
  cout << "Starting smoothing calculations..." << endl;
  // Smoothing step. This step updates probability distributions such that information in both past, present and future observations are incorporated
  smoo(nt) = phi(nt); // Last reconstruction is also a smoothed state
  for(int t=(nt-1); t>=1; --t){
    smoo(t) = elem_div(smoo(t+1),pred(t+1)+1.0e-12); 
    gridT = trans(smoo(t));
    for(int j=1; j<=ngx; ++j){
      gridT(j) = funvec(gridT(j),trans(Py(j)));
    }
    grid = trans(gridT);
    grid = funmat(grid,trans(Px));
    smoo(t) = value(grid);
    smoo(t) = elem_prod(smoo(t),phi(t));
  }
  report << smoo; 
  cout << "... finished smoothing calculations!" << endl;
  
TOP_OF_MAIN_SECTION
  // Allocate memory, required memory depends on grid resolution and length of data set
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1700000000); 
