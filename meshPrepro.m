function [meshGMF] = meshPrepro(meshGMF)
	% Victorien Menier Feb 2016
	% This function removes non-connected vertices
	% (useful if the mesh was generated using gmsh)
	
	% gmsh outputs 2D meshes as 3D meshes
	% -> Force dimension to be 2
	meshGMF.Dim = 2;
	
	% 1/ Tag the vertices that belong to a triangle
	% 2/ Compress vertices
	% 3/ Update vertex indices in the triangle and edge struct
	
	tag = zeros(1,meshGMF.NbrVer)';
	
	if meshGMF.Dim == 2
		%--- Tag vertices that belong to a trianle
		idx=reshape(meshGMF.Tri(:,1:3),[1,3*meshGMF.NbrTri]);
		tag(idx) = 1;
	end
	
	idx = find(tag==1);
	idx = sort(idx);
	
	NbrVerNew = size(idx);
	vidNew = [1:NbrVerNew]';
	
	tag(idx)=vidNew;
		
	%--- Compress vertices	
	
	meshGMF.Ver = meshGMF.Ver(idx,:);
	meshGMF.NbrVer = NbrVerNew;
	
	%--- Update triangles
	
	idx=reshape(meshGMF.Tri(:,1:3),[1,3*meshGMF.NbrTri]);
	idx=tag(idx);
	meshGMF.Tri(:,1:3) = reshape(idx,[meshGMF.NbrTri,3]);
	
	%--- Update Edges
	
	if ( meshGMF.NbrEfr > 0 )
		idx=reshape(meshGMF.Efr(:,1:2),[1,2*meshGMF.NbrEfr]);
		idx=tag(idx);
		meshGMF.Efr(:,1:2) = reshape(idx,[meshGMF.NbrEfr,2]);
	end
	
	%--- Remove edges of ref 11 (internal edges used for postprocessing)
	
	idx=find(meshGMF.Efr(:,3)~=11);
	meshGMF.Efr = meshGMF.Efr(idx,:);
	meshGMF.NbrEfr = size(idx,1);
	
end