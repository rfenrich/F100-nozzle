function [meshSU2] = convertGMFtoSU2 (meshGMF)
	
	% Victorien Menier Feb 2016
	% Convert mesh from the GMF (Inria) data structure to SU2
	% Useful, because boundary elements are dealt with differently for each forma
	
	meshSU2.dim   = meshGMF.Dim;
	meshSU2.npoin = meshGMF.NbrVer;
	meshSU2.poin  = [0];
	meshSU2.nelem = 0;
	meshSU2.elem  = [0];
	meshSU2.nmark = 0;
	meshSU2.mark(1).tag  = '';
	meshSU2.mark(1).nelem = 0;
	meshSU2.mark(1).elem  = [0];
	
	
	if meshGMF.Dim == 2
		meshSU2.nelem = meshGMF.NbrTri;
	else
		meshSU2.nelem = meshGMF.NbrTet;
	end
	
	
	% Vertices
	meshSU2.poin = [meshGMF.Ver(:,1:2), [0:1:meshGMF.NbrVer-1]'];
	
	% Tri, tet
	meshSU2.elem = [5*ones(meshGMF.NbrTri,1), meshGMF.Tri(:,1:3)-1, [0:1:meshGMF.NbrTri-1]'];
	
	% Bdry elt
	if meshGMF.NbrEfr > 0 
		ref = unique(meshGMF.Efr(:,3));
  	
		meshSU2.nmark = size(ref,1);
		
		for i=1:meshSU2.nmark
			meshSU2.mark(i).tag  = sprintf('%d', ref(i));
			idx = (meshGMF.Efr(:,3)==ref(i));
			meshSU2.mark(i).elem  = meshGMF.Efr(idx,1:2)-1;
			meshSU2.mark(i).nelem = size(meshSU2.mark(i).elem,1);
			meshSU2.mark(i).elem  = [3*ones(meshSU2.mark(i).nelem,1), meshSU2.mark(i).elem];
		end
	end
		
end