
function [flag] = checkCFDConvergence (ResFilNam)
	
	% checkCFDConvergence
	% Victorien Menier, Apr 2016
	% Open the residual history data file  (ResFilNam)
	% Check the final density residual R
	% If the final residual is higher than the initial one, the simulation is considered to have diverged
	% NB: If needed, more elaborate checks will be implemented in the future
	
	% Returns flag: false (0) if not converged, true (1) if converged	
	% If any error appears, true is returned.
	% -> false is only return if we know for sure that the computation has diverged
	
	
	% --- Open residual file and get header
	flag=0;
	dat=[];
	
	fid = fopen(ResFilNam,'r');
	
	if fid == -1
	  fprintf('  ## ERROR checkCFDConvergence : Open error. Check aborted, return true anyway.\n');
		flag=1;
		return
	end
	
	varNam=[];
	
	while ( ~feof(fid) )
		ll = fgetl(fid);
		
		varNam = strsplit(ll,{' ', ',', '"', '='});

		if (strcmp(varNam{1},'VARIABLES')) 
			break;
		end
	end
	
	if (feof(fid))
		fprintf('  ## ERROR checkCFDConvergence : Unexpected residual file header. Return true anyway.\n');
		flag=1;
		return
	end
	
	NbrVar = size(varNam,2);
	
	iRes=-1;
	for i=1:NbrVar
		if ( strcmp(varNam{i},'Res_Flow[0]') ) 
			iRes=i;
		end
	end
	
	if iRes <= 0
	  fprintf('  ## ERROR checkCFDConvergence : Unexpected residual file header (2). Return true anyway.\n');
		flag=1;
		return
	end
	
	ll = fgetl(fid); % Skip ZONE = ...
	ll = fgetl(fid);
		
	if (feof(fid) || strcmp(ll,'') == 1 )
		fprintf('  ## ERROR checkCFDConvergence : Unexpected residual file. Return true anyway.\n');
		flag=1;
		return
	end
	
	
	var = strsplit(ll,{' ', ',', '"', '='});
	resIni = str2double(var{iRes});
	
	while ( ~feof(fid) )
		ll = fgetl(fid);
	end
	
	var = strsplit(ll,{' ', ',', '"', '='});
	res = str2double(var{iRes});
	
	fprintf('  -- Check SU2 convergence.\n');
	fprintf('      Residual (density) = %f (ini) %f (final) ', resIni, res);
	if ( res > resIni ) 
		fprintf(' -> NOT converged.\n');
		flag=0;
	else
		fprintf(' -> convergence OK.\n');	
		flag=1;
	end

	fclose(fid);
	
end