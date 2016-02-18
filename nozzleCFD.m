function [ nozzle ] = nozzleCFD( fluid, freestream, nozzle, error, boundaryCdt )
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
	
	
	
	
	% ========================== GAS PROPERTIES ==============================
	gam = fluid.gam;
	R = fluid.R;
	
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
	
	% Calculate pressure ratio which determines state of nozzle:
	pressureRatio = inlet.Pstag/freestream.P;
	
	%% ======================= Get Nozzle Parameterization ==============
	
	%--- Inner wall
	[ A, dAdx, D ] = nozzleParameterization( nozzle );
	
	%--- Wall thickness
	t = @(x) piecewiseLinearGeometry(x,'t',nozzle.wall);
	
	%% ======================= MESH GENERATION ===========================
	
	fprintf('  -- Info : Mesh generation.\n');
	
	% Mesh size?	
	nozzle.sizWal = 0.01; % edge size around the nozzle wall
	nozzle.sizFar = 0.4;  % max size for the farfield region
	nozzle.sizSym = 0.1; % max size for the symmetry border
	
	xPosition = linspace(0,nozzle.geometry.length,100);
	fprintf('           axinoz.geo (input file to gmsh) created.\n');
	%nozzleCFDGmsh(nozzle, xPosition, A(xPosition') )
	nozzleCFDGmsh(nozzle, xPosition, D(xPosition')/2 )
	
	if(exist('axinoz.mesh', 'file') == 2)
	  delete('axinoz.mesh'); 
	end
	
	fprintf('           Calling gmsh (Cf gmsh.job).\n');
	!gmsh axinoz.geo -2 -o axinoz.mesh > gmsh.job
	
	if(exist('axinoz.mesh', 'file') ~= 2)
	  error('  ## ERROR : gmsh failed to generate the mesh. See gmsh.job for more details.'); 
	end
	
	
	fprintf('           %% %s created.\n', 'axinoz.mesh');
	
	%--- Convert to SU2 file format
	
	if(exist('axinoz.su2', 'file') == 2)
	  delete('axinoz.su2'); 
	end
	
	fprintf('           Mesh pre-processing/conversion to .su2\n');
	meshGMF = ReadGMF('axinoz.mesh');
	meshGMF = meshPrepro(meshGMF);
	meshSU2 = convertGMFtoSU2(meshGMF);
	WriteSU2(meshSU2, 'axinoz.su2');
	
	if(exist('axinoz.su2', 'file') ~= 2)
	  fprintf('  ## ERROR : MESH CONVERSION FROM .mesh to .su2 FAILED.\n');
		return;
	end
	
	fprintf('           %% %s created.\n', 'axinoz.su2');
	
	% ======================= RUN CFD SIMULATION (SU2) ===============
	
	% Write data file (.cfg) for SU2
	
	writeSU2DataFile( boundaryCdt );
	
	if(exist('history.dat', 'file') == 2)
	  delete('./history.dat');
	end
	
	if(exist('restart_flow.dat', 'file') == 2)
	  delete('./restart_flow.dat');
	end
	disp('  -- Info: running SU2 (Cf SU2.job )')
	!SU2_CFD euler_nozzle.cfg >SU2.job
	
	if(exist('restart_flow.dat', 'file') ~= 2)
	  disp('  ## ERROR nozzleCFD : SU2 simulation failed.')
	  return
	end
	
	% ======================= POST PROCESSING ===============
	
	% Extract averaged solution values at the nozzle's exit
	
	nozzle = nozzleCFDPostPro('restart_flow.dat', nozzle);
	
	massFlowRate = @(Pstag,Area,Tstag,M) (gam/((gam+1)/2)^((gam+1)/(2*(gam-1))))*Pstag*Area*AreaMachFunc(gam,M)/sqrt(gam*R*Tstag);
	nozzle.exit.Pstag = nozzle.inlet.Pstag*nozzle.PstagRatio;
	nozzle.massFlowRate = massFlowRate(nozzle.inlet.Pstag,nozzle.inlet.A,nozzle.inlet.Tstag,nozzle.flow.M(1));
	nozzle.approxThrust = nozzle.massFlowRate*(nozzle.exit.U - freestream.U) + (nozzle.exit.P - freestream.P)*nozzle.exit.A;
	
	fprintf('  CFD thrust = %f\n', nozzle.approxThrust);

end
