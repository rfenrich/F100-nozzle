function [] = solutionExtraction(Msh)
	
	ComputeR = @(v,v0,v1) ( ((v0(2)-v(2))*(v0(2)-v1(2)) - (v0(1)-v(1))*(v1(1)-v0(1) ) )/norm(v1-v0)^2);
	ComputeS = @(v,v0,v1) ( ((v0(2)-v(2))*(v1(1)-v0(1)) - (v0(1)-v(1))*(v1(2)-v0(2) ) )/norm(v1-v0)^2);
		
	v0=[0,0];
	v1=[0,2];
	
	v=[0.2,0.2];
	
	r=ComputeR(v,v0,v1)
	s=ComputeS(v,v0,v1)
	
	Msh.Ver(:,1:2)
	
	R=ComputeR(Msh.Ver(1,1:2),v0,v1)
	
	e2v=[1,2;2,3;3,1];
	
	%Tri=Msh.Tri;
	
	%Msh.Tri(:,e2v(1,:))
	
	%tmp=Msh.Ver(Msh.Tri(:,e2v(:,1)),:)
	%R = ComputeR(Msh.Ver(Msh.Tri(i,e2v(:,1)),:),v0,v1);
	
	R=[];
	S=[];
	
	for i=1:Msh.NbrTri
		for e=1:3
			R(e,1)=ComputeR(Msh.Ver(Msh.Tri(i,e2v(:,1)),:),v0,v1);
			S(e,1)=ComputeS(Msh.Ver(Msh.Tri(i,e2v(:,1)),:),v0,v1);
			R(e,2)=ComputeR(Msh.Ver(Msh.Tri(i,e2v(:,2)),:),v0,v1);
			S(e,2)=ComputeS(Msh.Ver(Msh.Tri(i,e2v(:,2)),:),v0,v1);
		end
		
	%	if ( )
		
		return
	end
	
end
