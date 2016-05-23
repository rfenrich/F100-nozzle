function [ nozzle ] = nozzleCFD( fluid, freestream, nozzle, err )
	% Victorien Menier Feb 2016
	% INPUTS:
	% nozzle geometry (i.e. border + wall thickness functions) 
	% Flow conditions:
	%  - Reference Mach number, Total Pressure Pt_ref, Total Temperature Tt_ref
	%  - Inlet : Total Temp Tt_in, Total Pres Pt_in
	%  - Outlet static pressure: Ps_out 
	% 
	% OUTPUTS:
	% nozzle = modified input structure with additional fields including flow
	% and specific geometry
	
	fprintf('\n\n--------------------------------------------------\n')
	fprintf('------------------  NOZZLECFD ------------------\n')
	fprintf('--------------------------------------------------\n\n')		
	
	% --- Hack used for the interface w/ dakota
	%      cf ./dakota-directory
	
	if ( exist('nozzle.workdir') )
		cd nozzle.workdir;
		addpath /Users/menier/calcul/DARPA/F100-nozzle;
	end
	
	% ---
	
	nozzle.output= 'PARAVIEW'; % TECPLOT OR PARAVIEW
	
	nozzle.success = 0;
	
	postpro = 0; % R
	
	% ========================== GAS PROPERTIES ==============================
	gam = fluid.gam;
	R   = fluid.R;
	
	% Area-Mach function from mass 1-D mass conservation equations:
	AreaMachFunc = @(g,M) ((g+1)/2)^((g+1)/(2*(g-1)))*M./(1+(g-1)*M.^2/2).^((g+1)/(2*(g-1)));
	% Sutherland's law for dynamic viscosity:
	dynamicViscosity = @(T) 1.716e-5*(T/273.15).^1.5*(273.15 + 110.4)./(T + 110.4); % kg/m*s
	
	% ========================= NOZZLE PROPERTIES ============================
	% Set inlet properties
	inlet = nozzle.inlet;
	
	% Calculate nozzle inlet, throat, and exit areas if they are not given:
	if(~exist('nozzle.inlet.A','var'))
	    nozzle.inlet.A = pi*inlet.D^2/4;
	end
	if(~exist('nozzle.throat.A','var'))
	    nozzle.throat.A = nozzle.inlet.A/nozzle.geometry.Ainlet2Athroat;
	end
	if(~exist('nozzle.exit.A','var'))
	    nozzle.exit.A = nozzle.geometry.Aexit2Athroat*nozzle.inlet.A/nozzle.geometry.Ainlet2Athroat;
	end
	
	nozzle.geometry.xApparentThroat = nozzle.geometry.xThroat; % initialize apparent throat location
	
	%% ======================= Get Nozzle Parameterization ==============
	
	%--- Inner wall
	[ A, dAdx, D, nozzle ] = nozzleParameterization( nozzle );
	
	%--- Wall thickness
	t = @(x) piecewiseLinearGeometry(x,'t',nozzle.wall);

	nozzle.geometry.DExit = D(nozzle.geometry.length);
    
    % --- Calc nozzle material volume
    % Volume calculation only works for spline parameterized nozzle geometry
    % and piecewise-linear parameterized nozzle wall thickness
    if(exist('pp','var')) % Exact volume for cubic spline parameterization
        nozzle.geometry.volume = wallVolume(pp,nozzle.wall);
    else % Approximate volume using trapezoidal integration
        nozzle.geometry.volume = wallVolume(nozzle.geometry.length,D,t,'integrate');
    end

	% ========================== CFD : boundary conditions
	
	nozzle.boundaryCdt.Mref  = freestream.M;
    nozzle.boundaryCdt.PsRef = freestream.P;
    nozzle.boundaryCdt.TsRef = freestream.T;
    nozzle.boundaryCdt.TtIn  = nozzle.inlet.Tstag;
    nozzle.boundaryCdt.PtIn  = nozzle.inlet.Pstag;
	nozzle.boundaryCdt.Uref  = freestream.U;
	nozzle.boundaryCdt.MuRef = dynamicViscosity(freestream.T);

	%--- RANS
	nozzle.boundaryCdt.LRey = D(nozzle.geometry.length)/2;
	
	nozzle.boundaryCdt.RhoRef = freestream.P/(R*freestream.T);
	nozzle.boundaryCdt.Re     = nozzle.boundaryCdt.RhoRef*freestream.U*nozzle.boundaryCdt.LRey/nozzle.boundaryCdt.MuRef;
	
	%% ======================= MESH GENERATION ===========================
	
	fprintf('	-- Mesh generation.\n');
	
	if ( strcmp(nozzle.meshSize,'coarse') ) 
		% Mesh size? -> Cf NozzleCFDGmsh()	
		nozzle.sizWal = 0.008;  % edge size around the nozzle wall
		nozzle.sizFar = 2;    % max size for the farfield region
		nozzle.sizSym = 0.08;   % max size for the symmetry border
		nozzle.yplus  = 2;      % y+ -> governs the minimal size of the 1st layer of the boundary layer mesh
	elseif ( strcmp(nozzle.meshSize,'medium') ) 
		nozzle.sizWal = 0.005;
		nozzle.sizFar = 1; 
		nozzle.sizSym = 0.05; 
		nozzle.yplus  = 1.5;
	elseif ( strcmp(nozzle.meshSize,'fine') ) 
		nozzle.sizWal = 0.0025; 
		nozzle.sizFar = 1; 
		nozzle.sizSym = 0.025; 
		nozzle.yplus  = 1;
	else 
		error ('  Unknown mesh complexity! (Must be either coarse, medium or fine)')
	end
	
	xPosition = linspace(0,nozzle.geometry.length,100);
	fprintf('		axinoz.geo (input file to gmsh) created.\n');
	%nozzleCFDGmsh(nozzle, xPosition, A(xPosition') )
	nozzleCFDGmsh(nozzle, xPosition, D(xPosition')/2 )
	
	if ( postpro == 1 )
		if(exist('axinoz.mesh', 'file') ~= 2)
		  error('NO axinoz.mesh was found!');
		end
	else
		if(exist('axinoz.mesh', 'file') == 2)
		  delete('axinoz.mesh'); 
		end
		
		fprintf('		Calling gmsh (Cf gmsh.job).\n');
		!gmsh axinoz.geo -2 -o axinoz.mesh > gmsh.job
		
		if(exist('axinoz.mesh', 'file') ~= 2)
		  error('  ## ERROR : gmsh failed to generate the mesh. See gmsh.job for more details.'); 
		end
		
		fprintf('		%% %s created.\n', 'axinoz.mesh');
		
	end
	
	%--- Convert to SU2 file format
	
	if(exist('axinoz.su2', 'file') == 2)
	  delete('axinoz.su2'); 
	end
	
	fprintf('		Mesh pre-processing/conversion to .su2\n');
	meshGMF = ReadGMF('axinoz.mesh');
	meshGMF = meshPrepro(meshGMF);
	meshSU2 = convertGMFtoSU2(meshGMF);
	WriteSU2(meshSU2, 'axinoz.su2');
	
	if(exist('axinoz.su2', 'file') ~= 2)
	  fprintf('  ## ERROR : MESH CONVERSION FROM .mesh to .su2 FAILED.\n');
		return;
	end
	
	fprintf('		%% %s created.\n', 'axinoz.su2');
	
	
	if ( postpro ~= 1 )
		% ======================= RUN CFD SIMULATION (SU2) ===============
		% Write data file (.cfg) for SU2
		
		if ( strcmp(nozzle.governing,'rans') )
			% No safe mode for RANS yet (input parameters are already safe for RANS)
			error('  ## ERROR nozzleCFD : Unable to converge the CFD solution.');
		end

		nozzle.CFDSafeMode = 0;
		
		writeSU2DataFile( nozzle );

		if(exist('history.dat', 'file') == 2)
		  delete('./history.dat');
		end
		
		if(exist('history.csv', 'file') == 2)
		  delete('./history.csv');
		end
		
		if(exist('restart_flow.dat', 'file') == 2)
		  delete('./restart_flow.dat');
		end
		
		%if ( strcmp(nozzle.governing,'rans') )
		%	tmp = input('  ## WARNING ! You might want to use an euler flow computation for now.\n Note: a viscous mesh was generated.\n Do you want to continue? (y/n) \n', 's');
		%	if ( ~strcmp(tmp,'y') )
		%		fprintf('STOP\n');
		%		return;
		%	else
		%		fprintf('-> Continue\n');
		%	end
		%end
		
		disp('	-- Running SU2 (Cf SU2.job )')
		!SU2_CFD axinoz.cfg >SU2.job
		
		% --- Restart simulation if not converged
		
		% Check the final residual of the simulation.
		% If converged, continue.
		% If diverged, restart the sim using a safer (but slower) set of input parameters.
		
		
		if ( strcmp(nozzle.output,'PARAVIEW') )
			ResFilNam='history.csv';
		else
			ResFilNam='history.dat';
		end
				
		if ( exist('restart_flow.dat', 'file') ~= 2 || checkCFDConvergence (ResFilNam) ~= 1 )
			
			if ( strcmp(nozzle.governing,'rans') )
				% No safe mode for RANS yet (input parameters are already safe for RANS)
				error('  ## Error nozzleCFD : Unable to converge the CFD solution.');
			end
			
			% The CFD solution is not converged : restart with the safe mode on
			nozzle.CFDSafeMode = 1;
			writeSU2DataFile( nozzle );
			disp('	-- Running SU2 for the 2nd time using safer parameters (Cf SU2_safe.job )')
			!SU2_CFD axinoz.cfg >SU2_safe.job
			if ( exist('restart_flow.dat', 'file') ~= 2 || checkCFDConvergence (ResFilNam) ~= 1 )
				error('  ## Error nozzleCFD : Unable to converge the CFD solution.');
			end
		end
	
	end
	
	% ======================= POST PROCESSING ===============
	
	SolSU2 = ReadSU2Sol('restart_flow.dat');
	
	% Extract averaged solution values at the nozzle's exit
	nozzle = nozzleCFDPostPro(meshSU2, SolSU2, nozzle, fluid, freestream);
	
	nozzle.success = 1;

end
