function [] = WriteGMFStruc (FilNam, xq, yq, vq )
	
	% Victorien Menier, Feb 2016
	% Writes a structured mesh defined by xq, yq
	% and a solution file, according to the values of vq
	
	[Ni,Nj] = size(xq)
	
	mesh.Dim   = 2;
	
	mesh.NbrVer = 0; % vertices
	mesh.Ver    = [0];
	
	mesh.NbrTri = 0; % triangles
	mesh.Tri    = [0];
	
	mesh.NbrEfr = 0; % boundary edges
	mesh.Efr    = [0];
	
	mesh.NbrTet = 0; % tetrahedra
	mesh.Tet    = [0];
	
	mesh.NbrVer = Ni*Nj;
	mesh.Ver = [reshape(xq,[mesh.NbrVer],1), reshape(yq,[mesh.NbrVer],1), zeros(mesh.NbrVer,1)]
	
	mesh.Tri = zeros((Nj-1)*(Ni-1)*2,4);
	mesh.NbrTri=0;
	for i=2:Ni
		for j=2:Nj
			ind = (i-2)*Nj+j-1;
			mesh.NbrTri = mesh.NbrTri+1;
			mesh.Tri(mesh.NbrTri,:) = [ind, ind+Nj, ind+Nj+1,1];
			mesh.NbrTri = mesh.NbrTri+1;
			mesh.Tri(mesh.NbrTri,:) = [ind, ind+1+Nj, ind+1,1];			
		end
	end 
	mesh.Tri(1:10,:)
	
	WriteGMF(FilNam,mesh,reshape(vq,[mesh.NbrVer],1))

end