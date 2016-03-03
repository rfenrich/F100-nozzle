function [SolSU2] = ReadSU2Sol(fname)
	
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
	 %tmp = sscanf(ll, '%f', NbrVar);
	 %dat(i,:) = tmp';
     tmp = str2double(strsplit(ll));
     dat(i,:) = tmp(1:NbrVar);
	 i=i+1;
	end
	
	NbrVer=i-1;
	
	SolSU2.NbrVer = NbrVer;
	SolSU2.Dat    = dat;
	SolSU2.NamVar = varNam;
	
end