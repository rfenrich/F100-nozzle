function mesh = ReadGMF(name)
	% Victorien Menier Feb 2016
	% Read Inria Gamma Mesh Format (GMF)
	% For more information about GMF,
	% see https://www.rocq.inria.fr/gamma/gamma/Membres/CIPD/Loic.Marechal/Research/LM6.html
	% mesh = ReadGMF(meshname)
	% reads meshname.mesh into struct mesh
	
	%% --- Initialize

	mesh.Dim   = 0;

	mesh.NbrVer = 0; % vertices
	mesh.Ver    = [0];

	mesh.NbrTri = 0; % triangles
	mesh.Tri    = [0];

	mesh.NbrEfr = 0; % boundary edges
	mesh.Efr    = [0];

	mesh.NbrTet = 0; % tetrahedra
	mesh.Tet    = [0];
	
	
	fid = fopen(name,'r');
	if ( fid == -1 ) 
	 error(['Mesh file ' name ' does not exist']);
	%else
	% disp([ '% ' name ' OPENED ']);
	end
	
	
	while ( ~feof(fid) )
	
	  str = fgetl( fid );
	  str = str(find(str~=' '));
	  switch ( lower(str) )
	  
	  case 'dimension'
	    mesh.Dim = fscanf( fid, '%d', 1 );
	    if ( mesh.Dim ~= 2 & mesh.Dim ~=3 )
	      error(' Invalid inout mesh ');
	    end  
	
	  case 'vertices'
	    mesh.NbrVer = fscanf( fid, '%d', 1 );
	    if (  mesh.Dim  == 2 )
	      mesh.Ver = fscanf( fid, '%f %f %d', [3,mesh.NbrVer] );
	    else
	      mesh.Ver = fscanf( fid, '%f %f %f %d', [4,mesh.NbrVer] );
	    end
			mesh.Ver = mesh.Ver';
	   
	  case 'triangles'
	    mesh.NbrTri = fscanf( fid, '%d', 1 );
	    mesh.Tri = fscanf( fid, '%d %d %d %d', [4,mesh.NbrTri] );
			mesh.Tri = mesh.Tri';
	  case 'tetrahedra'
	    mesh.NbrTet = fscanf( fid, '%d', 1 );
			mesh.Tet = fscanf( fid, '%d %d %d %d %d', [5,mesh.NbrTri] );
	  	mesh.Tet = mesh.Tet'
	  case 'edges'
	    mesh.NbrEfr = fscanf( fid, '%d', 1 );
	    mesh.Efr = fscanf( fid, '%d %d %d', [3,mesh.NbrEfr] );  
			mesh.Efr = mesh.Efr';
	  end
	end
	
	
	fclose(fid);
end

