
function [nozzle] = nozzleCFDPostPro(SolNam, nozzle)
	
	[NbrVer,dat,namVar] = ReadSU2Sol(SolNam);
	
	NbrVar = size(namVar,2);
	
	xextract = 0.67; % x coordinate of the exit
	
	ix     = -1;
	iMach  = -1;
	iTem   = -1;
	iRho   = -1;
	iCons1 = -1;
	iCons2 = -1;
	iCons3 = -1;
	iCons4 = -1;
	iPres  = -1;
	
	for i=1:NbrVar
		if ( strcmp(namVar(i),'"x"') ) 
			ix = i;
		elseif ( strcmp(namVar(i),'"y"') ) 
			iy = i;
		elseif ( strcmp(namVar(i),'"Mach"') ) 
			iMach = i;
		elseif ( strcmp(namVar(i),'"Temperature"') ) 
			iTem = i;
		elseif ( strcmp(namVar(i),'"Conservative_1"') ) 
			iCons1 = i;
		elseif ( strcmp(namVar(i),'"Conservative_2"') ) 
			iCons2 = i;
		elseif ( strcmp(namVar(i),'"Conservative_3"') ) 
			iCons3 = i;
		elseif ( strcmp(namVar(i),'"Conservative_4"') ) 
			iCons4 = i;
		elseif ( strcmp(namVar(i),'"Pressure"') ) 
			iPres = i;
		end
	end
	
	idx = find( dat(:,ix) == xextract & dat(:,iy) < 0.2919 );
	DatLin = dat(idx,:);
	DatLin = sortrows(DatLin,iy);
	
	NbvLin = size(idx,1);
	
	if ( NbvLin < 2 ) 
		error('  ## ERROR nozzleCFDPostPro : No extraction data was found');
	end 
	
	avgDat=zeros(1,size(DatLin,2));
  
  for i=2:NbvLin
  	dy = DatLin(i,iy)-DatLin(i-1,iy);
  	avgDat(1,:) = avgDat(1,:)+dy*DatLin(i,:);
  end
	
	hei = DatLin(NbvLin,iy)-DatLin(1,iy);
	avgDat(1,:) = avgDat(1,:)/hei;
	
	%nozzle.exit.Pstag = nozzle.inlet.Pstag*nozzle.PstagRatio;
	nozzle.exit.M = avgDat(1,iMach);
	nozzle.exit.T = avgDat(1,iTem);
	nozzle.exit.U = (avgDat(1,iCons2)+avgDat(1,iCons3))/avgDat(1,iCons1);
	nozzle.exit.P = avgDat(1,iPres);
	
	fprintf('\n -- Info: CFD results (averaged values at nozzle exit)\n');
	fprintf('           Mach = %f\n', nozzle.exit.M );
	fprintf('           Temp = %f\n', nozzle.exit.T );
	fprintf('           Vel  = %f\n', nozzle.exit.U );
	fprintf('           Pres = %f\n', nozzle.exit.P );
	
end