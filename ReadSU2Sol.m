function [NbrVer,dat,varNam] = ReadSU2Sol(fname)
	
	dat=[];
	
	fid = fopen(fname,'r');
	
	if fid == -1
	  error('  ## ERROR ReadSU2Sol : Open error\n');
	end
	
	ll = fgetl(fid);
	varNam = strsplit(ll,'	');
	
	NbrVar = size(varNam,2);
	
	if ( NbrVar < 1 )
		error('  ## ERROR ReadSU2Sol : Invalid solution file.\n');
	end
	 
	i=1;
	while ( ~feof(fid) )
	 ll = fgetl(fid);
	 tmp = sscanf(ll, '%f', NbrVar);
	 dat(i,:) = tmp';
	 i=i+1;
	end
	
	NbrVer=i-1;
	
end