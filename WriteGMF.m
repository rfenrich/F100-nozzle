function ok = WriteGMF(FilNam,Msh,Sol)
	
	% Victorien Menier Feb 2016
	% Writes FilNam.mesh 
	% and FilNam.sol is Sol is provided
	% GMF : Gamma Mesh Format > Cf https://www.rocq.inria.fr/gamma/gamma/Membres/CIPD/Loic.Marechal/Research/LM6.html
	
	
	% --- Write mesh file
	
	k = strfind(FilNam,'.mesh'); 
  if ( length(k) ~= 0 )
    in  = [ FilNam(1:(k(1)-1)) '.mesh'];
  else
    in  = [ FilNam '.mesh'];
  end
	
  fid = fopen(FilNam,'w+');
	if ( fid == -1 ) 
    error(['Cannot create mesh file ' FilNam ]);
  else
    disp([ '% '  FilNam ' CREATED ']);
  end

	fprintf(fid,'MeshVersionFormatted 2\n');
	fprintf(fid,'Dimension\n%d\n', Msh.Dim);
	
	fprintf(fid,'Vertices\n%d\n',Msh.NbrVer);
	if Msh.Dim == 2 
		fprintf(fid,' %f %f %d \n', Msh.Ver');
	else
		fprintf(fid,' %f %f %f %d \n', Msh.Ver');
	end
	
	fprintf(fid,'Triangles\n%d\n',Msh.NbrTri);
	fprintf(fid,' %d %d %d  %d\n',Msh.Tri');
	fprintf(fid,'Edges\n%d\n',Msh.NbrEfr);
	fprintf(fid,' %d %d %d \n',Msh.Efr');
	
	fprintf(fid,'\nEnd\n');
	
	fclose(fid);
	
	% --- Write solution if one is provided
	
	[NbvSol, SolSiz] = size(Sol);
	
	if ( NbvSol ~= Msh.NbrVer )
		fprintf('  ## ERROR WriteGMF : INCORRECT SOLUTION FILE\n');
		return;
	end
	
	k = findstr(FilNam,'.mesh')
	
	if ( k < 2  )
		fprintf('  ## ERROR WriteGMF : INCORRECT FILE NAME\n');
		return;
	end	
	
	SolNam = [FilNam(1:k-1),'.sol'];
	
	fid = fopen(SolNam,'w+');
	if ( fid == -1 ) 
    error(['Cannot create solution file ' SolNam ]);
  else
    disp([ '% '  SolNam ' CREATED ']);
  end
	
	
	fprintf(fid,'MeshVersionFormatted 2\n');
	fprintf(fid,'Dimension\n%d\n', Msh.Dim);
	fprintf(fid,'SolAtVertices\n%d\n',Msh.NbrVer);

	
	fprintf(fid,'%d ',SolSiz);
	for i=1:SolSiz
		fprintf(fid,'1 ');
	end
	fprintf(fid,'\n');
	
	for i=1:Msh.NbrVer
		for j=1:SolSiz
			fprintf(fid,'%f ', Sol(i,j));
		end
		fprintf(fid,'\n');
	end
	
	fprintf(fid,'\nEnd\n');
	fclose(fid);
	

end