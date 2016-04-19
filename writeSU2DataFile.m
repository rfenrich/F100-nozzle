% 01-20-2015
% Victorien Menier

function [] = writeSU2DataFile ( nozzle )
  
  % TODO later:
  % Adapt the parameters to the stiffness of the problem (CFL, nbr of sub-iterations of the Newton method, etc.)

	rans = strcmp(nozzle.governing,'rans');
	
	safeMode =  nozzle.CFDSafeMode;
	
	if ( rans ) 
		% --- We use multigrid for RANS (not for Euler)
		%      -> fewer iterations required
		nozzle.NbrIte = 200;
	else if ( safeMode )
		nozzle.NbrIte = 600;
	else 
		nozzle.NbrIte = 300;
	end
	
  fprintf('	-- Writing SU2 datafile axinoz.cfg\n')
  fprintf('		Freestream Mach:               %f\n'   ,nozzle.boundaryCdt.Mref );
  fprintf('		Freestream Static Pressure:    %f Pa\n',nozzle.boundaryCdt.PsRef);
  fprintf('		Freestream Static Temperature: %f K\n' ,nozzle.boundaryCdt.TsRef);
  fprintf('		Inlet Stagnation Pressure:     %f Pa\n',nozzle.boundaryCdt.TtIn );
  fprintf('		Inlet Stagnation Temperature:  %f K\n' ,nozzle.boundaryCdt.PtIn );
  %fprintf('           Outlet Static Pressure:        %f Pa\n',nozzle.boundaryCdt.PsOut);
  
  DatOut=fopen('axinoz.cfg','w');
  
  fprintf(DatOut,'%% ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------%%  \n');
  fprintf(DatOut,'%% Physical governing equations (EULER, NAVIER_STOKES)                             \n');
  
	if ( ~rans )
		fprintf(DatOut,' PHYSICAL_PROBLEM= EULER                                                           \n');
  else
		fprintf(DatOut,' PHYSICAL_PROBLEM= NAVIER_STOKES                                                   \n');
		fprintf(DatOut,' KIND_TURB_MODEL= SA                                                               \n');
	end
	
	fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Mathematical problem (DIRECT, CONTINUOUS_ADJOINT)                               \n');
  fprintf(DatOut,' MATH_PROBLEM= DIRECT                                                              \n');
  fprintf(DatOut,'%% Restart solution (NO, YES)                                                      \n');
  fprintf(DatOut,' RESTART_SOL= NO                                                                   \n');
  fprintf(DatOut,' SYSTEM_MEASUREMENTS= SI                                                           \n');
  fprintf(DatOut,'%% -------------------- COMPRESSIBLE FREE-STREAM DEFINITION --------------------%% \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Mach number (non-dimensional, based on the free-stream values)                  \n');
  fprintf(DatOut,' MACH_NUMBER= %f                                                                   \n', nozzle.boundaryCdt.Mref);
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Angle of attack (degrees)                                                       \n');
  fprintf(DatOut,' AoA= 0.0                                                                          \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Side-slip angle (degrees)                                                       \n');
  fprintf(DatOut,' SIDESLIP_ANGLE= 0.0                                                               \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Free-stream pressure (101325.0 N/m^2 by default, only for Euler equations)      \n');
  fprintf(DatOut,'FREESTREAM_PRESSURE= %f                                                            \n', nozzle.boundaryCdt.PsRef);
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Free-stream temperature (288.15 K by default)                                   \n');
  fprintf(DatOut,'FREESTREAM_TEMPERATURE= %f                                                         \n', nozzle.boundaryCdt.TsRef);

	if ( rans )                                                                                        
		fprintf(DatOut,'REYNOLDS_NUMBER= 2154263.930611                                                  \n');
		fprintf(DatOut,'%                                                                                \n');
		fprintf(DatOut,'% Reynolds length (1 m by default)                                               \n');
		fprintf(DatOut,'REYNOLDS_LENGTH= 0.304800                                                        \n');
	end

  fprintf(DatOut,'                                                                                   \n');
  fprintf(DatOut,'%% ---------------------- REFERENCE VALUE DEFINITION ---------------------------%%  \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Reference origin for moment computation                                         \n');
  fprintf(DatOut,' REF_ORIGIN_MOMENT_X = 0.25                                                        \n');
  fprintf(DatOut,' REF_ORIGIN_MOMENT_Y = 0.00                                                        \n');
  fprintf(DatOut,' REF_ORIGIN_MOMENT_Z = 0.00                                                        \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Reference length for pitching, rolling, and yaMAIN_BOX non-dimensional moment   \n');
  fprintf(DatOut,' REF_LENGTH_MOMENT= 1.0                                                            \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Reference area for force coefficients (0 implies automatic calculation)         \n');
  fprintf(DatOut,' REF_AREA= 0                                                                       \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Flow non-dimensionalization (DIMENSIONAL, FREESTREAM_PRESS_EQ_ONE,              \n');
  fprintf(DatOut,'%%                              FREESTREAM_VEL_EQ_MACH, FREESTREAM_VEL_EQ_ONE)     \n');
  fprintf(DatOut,' REF_DIMENSIONALIZATION= DIMENSIONAL                                     \n');
  fprintf(DatOut,'                                                                                    \n');
  fprintf(DatOut,'%% ----------------------- BOUNDARY CONDITION DEFINITION -----------------------%%  \n');
  fprintf(DatOut,'%%                                                                                  \n');
  fprintf(DatOut,'%% Marker of the Euler boundary (0 implies no marker)                               \n');

	if ( ~rans )
  	fprintf(DatOut,' MARKER_EULER= ( 6, 7, 8, 9 )                                                     \n');
	else 
		fprintf(DatOut,'  MARKER_HEATFLUX= ( 6, 0.0, 7, 0.0, 8, 0.0, 9 , 0.0 )                            \n');
	end
	
	fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'                                                                                   \n');
  fprintf(DatOut,'%% Inlet boundary marker(s) (NONE = no marker)                                     \n');
  fprintf(DatOut,'%% Format: ( inlet marker, total temperature, total pressure, flow_direction_x,    \n');
  fprintf(DatOut,'%%           flow_direction_y, flow_direction_z, ... )                             \n');

	% INLET CONDITIONS
	
  fprintf(DatOut,' MARKER_INLET= ( 10, %f, %f, 1.0, 0.0, 0.0 ) \n ', nozzle.boundaryCdt.TtIn, nozzle.boundaryCdt.PtIn);

	%fprintf(DatOut,'  4, %f, %f, 1.0, 0.0, 0.0,  ', nozzle.boundaryCdt.TsRef, nozzle.boundaryCdt.PsRef);
	%fprintf(DatOut,'  5, %f, %f, 1.0, 0.0, 0.0 ) \n', nozzle.boundaryCdt.TsRef, nozzle.boundaryCdt.PsRef);
  fprintf(DatOut,' MARKER_FAR= ( 4, 5 )                                                               \n');
	
  fprintf(DatOut,'%% Marker of symmetry boundary (0 implies no marker)                                \n');
  fprintf(DatOut,' MARKER_SYM= ( 1, 2 )                                                               \n');
  fprintf(DatOut,'                                                                                    \n');
  fprintf(DatOut,' MARKER_OUTLET= ( 3, %f)                                                            \n', nozzle.boundaryCdt.PsRef);
  fprintf(DatOut,'%% = 0.89 * Ptot                                                                    \n');
  fprintf(DatOut,'                                                                                    \n');
  fprintf(DatOut,'%% ------------- COMMON PARAMETERS TO DEFINE THE NUMERICAL METHOD --------------%%  \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Numerical method for spatial gradients (GREEN_GAUSS, WEIGHTED_LEAST_SQUARES)    \n');
  fprintf(DatOut,' NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES                                           \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Objective function in optimization problem (DRAG, LIFT, SIDEFORCE, MOMENT_X,    \n');
  fprintf(DatOut,'%%                                             MOMENT_Y, MOMENT_Z, EFFICIENCY,     \n');
  fprintf(DatOut,'%%                                             EQUIVALENT_AREA, NEARFIELD_PRESSURE,\n');
  fprintf(DatOut,'%%                                             FORCE_X, FORCE_Y, FORCE_Z, THRUST,  \n');
  fprintf(DatOut,'%%                                             TORQUE, FREE_SURFACE, TOTAL_HEATFLUX\n');
  fprintf(DatOut,'%%                                             MAXIMUM_HEATFLUX, INVERSE_DESIGN_PRE\n');
  fprintf(DatOut,'%%                                             INVERSE_DESIGN_HEATFLUX)            \n');
  fprintf(DatOut,'%% OBJECTIVE_FUNCTION= DRAG                                                          \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Courant-Friedrichs-Lewy condition of the finest grid                            \n');
	if ( safeMode )
  	fprintf(DatOut,' CFL_NUMBER= 25                                                                  \n');
	else
		fprintf(DatOut,' CFL_NUMBER= 25                                                                  \n');
	end
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Adaptive CFL number (NO, YES)                                                   \n');
  fprintf(DatOut,' CFL_ADAPT= NO                                                                    \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Parameters of the adaptive CFL number (factor down, factor up, CFL min value,   \n');
  fprintf(DatOut,'%%                                        CFL max value )                          \n');
  fprintf(DatOut,' CFL_ADAPT_PARAM= ( 1.5, 0.5, 1.0, 100.0 )                                         \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Runge-Kutta alpha coefficients                                                  \n');
  fprintf(DatOut,' RK_ALPHA_COEFF= ( 0.66667, 0.66667, 1.000000 )                                    \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Number of total iterations                                                      \n');
  fprintf(DatOut,' EXT_ITER= %d                                                                 \n', nozzle.NbrIte );
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Linear solver for the implicit formulation (BCGSTAB, FGMRES)                    \n');
  fprintf(DatOut,' LINEAR_SOLVER= FGMRES                                                             \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Min error of the linear solver for the implicit formulation                     \n');
  fprintf(DatOut,' LINEAR_SOLVER_ERROR= 1E-6                                                         \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Max number of iterations of the linear solver for the implicit formulation      \n');
  fprintf(DatOut,' LINEAR_SOLVER_ITER= 20                                                             \n');
  fprintf(DatOut,'                                                                                    \n');
  fprintf(DatOut,'%% ----------------------- SLOPE LIMITER DEFINITION ----------------------------%%  \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Reference element length for computing the slope and sharp edges limiters.      \n');
  fprintf(DatOut,' REF_ELEM_LENGTH= 0.1                                                              \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Coefficient for the limiter                                                     \n');
  fprintf(DatOut,' LIMITER_COEFF= 0.3                                                                \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Coefficient for the sharp edges limiter                                         \n');
  fprintf(DatOut,' SHARP_EDGES_COEFF= 3.0                                                            \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Reference coefficient (sensitivity) for detecting sharp edges.                  \n');
  fprintf(DatOut,' REF_SHARP_EDGES= 3.0                                                              \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Remove sharp edges from the sensitivity evaluation (NO, YES)                    \n');
  fprintf(DatOut,' SENS_REMOVE_SHARP= YES                                                            \n');
  fprintf(DatOut,'                                                                                    \n');

	if ( rans )
  	fprintf(DatOut,'%% -------------------------- MULTIGRID PARAMETERS -----------------------------%%  \n');
  	fprintf(DatOut,'%%                                                                                 \n');
  	fprintf(DatOut,'%% Multi-Grid Levels (0 = no multi-grid)                                           \n');
  	fprintf(DatOut,' MGLEVEL= 3                                                                        \n');
  	fprintf(DatOut,'%%                                                                                 \n');
  	fprintf(DatOut,'%% Multi-grid cycle (V_CYCLE, W_CYCLE, FULLMG_CYCLE)                               \n');
  	fprintf(DatOut,' MGCYCLE= V_CYCLE                                                                  \n');
  	fprintf(DatOut,'%%                                                                                 \n');
  	fprintf(DatOut,'%% Multi-Grid PreSmoothing Level                                                   \n');
  	fprintf(DatOut,' MG_PRE_SMOOTH= ( 1, 2, 3, 3 )                                                     \n');
  	fprintf(DatOut,'%%                                                                                 \n');
  	fprintf(DatOut,'%% Multi-Grid PostSmoothing Level                                                  \n');
  	fprintf(DatOut,' MG_POST_SMOOTH= ( 0, 0, 0, 0 )                                                    \n');
  	fprintf(DatOut,'%%                                                                                 \n');
  	fprintf(DatOut,'%% Jacobi implicit smoothing of the correction                                     \n');
  	fprintf(DatOut,' MG_CORRECTION_SMOOTH= ( 0, 0, 0, 0 )                                              \n');
  	fprintf(DatOut,'%%                                                                                 \n');
  	fprintf(DatOut,'%% Damping factor for the residual restriction                                     \n');
  	fprintf(DatOut,' MG_DAMP_RESTRICTION= 0.9                                                          \n');
  	fprintf(DatOut,'%%                                                                                 \n');
  	fprintf(DatOut,'%% Damping factor for the correction prolongation                                  \n');
  	fprintf(DatOut,' MG_DAMP_PROLONGATION= 0.9                                                         \n');
  	fprintf(DatOut,'                                                                                   \n');
	end
	
  fprintf(DatOut,'%% -------------------- FLOW NUMERICAL METHOD DEFINITION -----------------------%% \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Convective numerical method (JST, LAX-FRIEDRICH, CUSP, ROE, AUSM, HLLC,         \n');
  fprintf(DatOut,'%%                              TURKEL_PREC, MSW)                                  \n');
  fprintf(DatOut,' CONV_NUM_METHOD_FLOW= ROE                                                         \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Spatial numerical order integration (1ST_ORDER, 2ND_ORDER, 2ND_ORDER_LIMITER)   \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,' SPATIAL_ORDER_FLOW= 2ND_ORDER_LIMITER                                              \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Slope limiter (VENKATAKRISHNAN, MINMOD)                                         \n');
  fprintf(DatOut,' SLOPE_LIMITER_FLOW= VENKATAKRISHNAN                                               \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% 1st, 2nd and 4th order artificial dissipation coefficients                      \n');
  fprintf(DatOut,' AD_COEFF_FLOW= ( 0.15, 0.5, 0.02 )                                                \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Time discretization (RUNGE-KUTTA_EXPLICIT, EULER_IMPLICIT, EULER_EXPLICIT)      \n');
  fprintf(DatOut,' TIME_DISCRE_FLOW= EULER_IMPLICIT                                                  \n');
  fprintf(DatOut,'                                                                                   \n');
	
	if ( rans ) 
		
		fprintf(DatOut,'% -------------------- TURBULENT NUMERICAL METHOD DEFINITION ------------------%   \n');
		fprintf(DatOut,'% Convective numerical method (SCALAR_UPWIND)                                      \n');
		fprintf(DatOut,'CONV_NUM_METHOD_TURB= SCALAR_UPWIND                                                \n');                                                                               
		fprintf(DatOut,'% Spatial numerical order integration (1ST_ORDER, 2ND_ORDER, 2ND_ORDER_LIMITER)    \n');
		fprintf(DatOut,'SPATIAL_ORDER_TURB= 1ST_ORDER                                                      \n');
		fprintf(DatOut,'TIME_DISCRE_TURB= EULER_IMPLICIT                                                   \n');	
	end
	
	if ( rans || safeMode ) 
		fprintf(DatOut,'RELAXATION_FACTOR_FLOW= 0.1  \n');
	end
	
  fprintf(DatOut,'%% --------------------------- CONVERGENCE PARAMETERS --------------------------&  \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Convergence criteria (CAUCHY, RESIDUAL)                                         \n');
  fprintf(DatOut,' CONV_CRITERIA= RESIDUAL                                                           \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Residual reduction (order of magnitude with respect to the initial value)       \n');
  fprintf(DatOut,' RESIDUAL_REDUCTION= 8                                                             \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Min value of the residual (log10 of the residual)                               \n');
  fprintf(DatOut,' RESIDUAL_MINVAL= -12                                                              \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Start convergence criteria at iteration number                                  \n');
  fprintf(DatOut,' STARTCONV_ITER= 25                                                                \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Number of elements to apply the criteria                                        \n');
  fprintf(DatOut,' CAUCHY_ELEMS= 100                                                                 \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Epsilon to control the series convergence                                       \n');
  fprintf(DatOut,' CAUCHY_EPS= 1E-10                                                                 \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Function to apply the criteria (LIFT, DRAG, NEARFIELD_PRESS, SENS_GEOMETRY,     \n');
  fprintf(DatOut,'%% 	      	    		 SENS_MACH, DELTA_LIFT, DELTA_DRAG)                           \n ');
  fprintf(DatOut,' CAUCHY_FUNC_FLOW= DRAG                                                            \n');
  fprintf(DatOut,'                                                                                    \n');
  fprintf(DatOut,'%% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%%  \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Mesh input file                                                                 \n');
  fprintf(DatOut,' MESH_FILENAME= axinoz.su2                                                         \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Mesh output file                                                                \n');
  fprintf(DatOut,' MESH_OUT_FILENAME= mesh_out.su2                                                   \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Restart flow input file                                                         \n');
  fprintf(DatOut,' SOLUTION_FLOW_FILENAME= solution_flow.dat                                         \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Restart adjoint input file                                                      \n');
  fprintf(DatOut,' SOLUTION_ADJ_FILENAME= solution_adj.dat                                           \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Mesh input file format (SU2)                                                    \n');
  fprintf(DatOut,' MESH_FORMAT= SU2                                                                  \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Output file format (PARAVIEW, TECPLOT)                                          \n');
  fprintf(DatOut,' OUTPUT_FORMAT= TECPLOT                                                            \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Output file convergence history                                                 \n');
  fprintf(DatOut,' CONV_FILENAME= history                                                            \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Output file restart flow                                                        \n');
  fprintf(DatOut,' RESTART_FLOW_FILENAME= restart_flow.dat                                           \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Output file restart adjoint                                                     \n');
  fprintf(DatOut,' RESTART_ADJ_FILENAME= restart_adj.dat                                             \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Output file flow (w/o extension) variables                                      \n');
  fprintf(DatOut,' VOLUME_FLOW_FILENAME= flow                                                        \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Output file adjoint (w/o extension) variables                                   \n');
  fprintf(DatOut,' VOLUME_ADJ_FILENAME= adjoint                                                      \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Output Objective function gradient (using continuous adjoint)                   \n');
  fprintf(DatOut,' GRAD_OBJFUNC_FILENAME= of_grad.dat                                                \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Output file surface flow coefficient (w/o extension)                            \n');
  fprintf(DatOut,' SURFACE_FLOW_FILENAME= surface_flow                                               \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Output file surface adjoint coefficient (w/o extension)                         \n');
  fprintf(DatOut,' SURFACE_ADJ_FILENAME= surface_adjoint                                             \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Writing solution frequency                                                      \n');
  fprintf(DatOut,' WRT_SOL_FREQ= 1000                                                                 \n');
  fprintf(DatOut,'%%                                                                                 \n');
  fprintf(DatOut,'%% Writing convergence history frequency                                           \n');
  fprintf(DatOut,' WRT_CON_FREQ= 1                                                                   \n');

  fclose(DatOut);
end