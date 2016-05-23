function [ A, dAdx, D, nozzle ] = nozzleParameterization( nozzle )
	if(strcmp(nozzle.geometry.shape,'spline'))
	    % set up spline for nozzle
	    if(ischar(nozzle.geometry.spline.seed)) % seed shape is given
	        % Set control points for splines (xNode), value of function at control
	        % point (yNode), and slopes at start and end of spline (slopes)
	        xNode = nozzle.geometry.spline.breaks;
	        if(max(xNode) > nozzle.geometry.length || min(xNode) < 0) % check user given location of control points
	            error('Spline control point outside nozzle length domain');
	        end
	        yNode = nozzleGeometry(xNode,'D',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.spline.seed)/2;
	        nozzle.geometry.spline.seed = [xNode, yNode];
	    else
	        % Extract control points for splines and their values from the
	        % given array
	        xNode = nozzle.geometry.spline.seed(:,1);
	        yNode = nozzle.geometry.spline.seed(:,2);
	        if(xNode(1) == 0 && xNode(end) == nozzle.geometry.length) % check user given control point values (yNode) match user given area ratios
	            areaRatio = yNode(end)^2/yNode(1)^2;
	            areaRatioTolerance = 1e-3;
	            if(areaRatio > nozzle.geometry.Aexit2Athroat/nozzle.geometry.Ainlet2Athroat + areaRatioTolerance || areaRatio < nozzle.geometry.Aexit2Athroat/nozzle.geometry.Ainlet2Athroat - areaRatioTolerance)
	               error('Spline control point values do not match given nozzle area ratios'); 
	            end
	        end
	    end
			
	    slopes = nozzle.geometry.spline.slopes;
	    pp = spline(xNode,[slopes(1); yNode; slopes(2)]); % perform piecewise cubic spline interpolation

	    % Adjust nozzle throat size/location information if it has changed
	    [xThroat, yThroat] = splineGeometry(0, 'throat', pp);
	    if(abs(xThroat - nozzle.geometry.xThroat)/xThroat >= 1e-6)
	        fprintf('throat size/location changed with spline parameterization\n');
	    end
	    nozzle.geometry.xThroat = xThroat;
	    nozzle.geometry.xApparentThroat = nozzle.geometry.xThroat; % initialize apparent throat location
	    nozzle.throat.A = pi*yThroat^2;
	    nozzle.geometry.Ainlet2Athroat = nozzle.inlet.A/nozzle.throat.A;
	    nozzle.geometry.Aexit2Athroat = nozzle.exit.A/nozzle.throat.A;

	    % Make necessary functions for splined nozzle shape
	    A = @(x) splineGeometry(x,'A',pp);
	    dAdx = @(x) splineGeometry(x,'dAdx',pp);
	    D = @(x) splineGeometry(x,'D',pp);

	elseif(strcmp(nozzle.geometry.shape,'B-spline'))

	    % Adjust nozzle throat size/location information if it has changed
	    if(nozzle.geometry.bSpline.degree == 2)
	        [xThroat, yThroat] = BsplineGeometry(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);
	    elseif(nozzle.geometry.bSpline.degree == 3)
	        [xThroat, yThroat] = BsplineGeometry3(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);
	    else
	        error('Only B-splines of degree 2 and 3 are supported');
	    end

	    if(xThroat ~= nozzle.geometry.xThroat)
	        fprintf('throat size/location changed with spline parameterization\n');
	    end
	    nozzle.geometry.xThroat = xThroat;
	    nozzle.geometry.xApparentThroat = nozzle.geometry.xThroat; % initialize apparent throat location
	    nozzle.throat.A = pi*yThroat^2;
        nozzle.inlet.A = pi*nozzle.geometry.bSpline.coefs(2,1)^2;
        nozzle.exit.A = pi*nozzle.geometry.bSpline.coefs(2,end)^2;        
	    nozzle.geometry.Ainlet2Athroat = nozzle.inlet.A/nozzle.throat.A;
	    nozzle.geometry.Aexit2Athroat = nozzle.exit.A/nozzle.throat.A;

	    % Make necessary functions for splined nozzle shape
	    if(nozzle.geometry.bSpline.degree == 2)
	        A = @(x) BsplineGeometry(x,'A',nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs);
	        dAdx = @(x) BsplineGeometry(x,'dAdx',nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs);
	        D = @(x) BsplineGeometry(x,'D',nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs);
	    elseif(nozzle.geometry.bSpline.degree == 3)
	        A = @(x) BsplineGeometry3(x,'A',nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs);
	        dAdx = @(x) BsplineGeometry3(x,'dAdx',nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs);
	        D = @(x) BsplineGeometry3(x,'D',nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs);
	    end

	elseif(strcmp(nozzle.geometry.shape,'B-spline-mex'))

	    % Adjust nozzle throat size/location information if it has changed
	    if(nozzle.geometry.bSpline.degree == 2)
	        [xThroat, yThroat] = BsplineGeometry(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);
	    elseif(nozzle.geometry.bSpline.degree == 3)
	        [xThroat, yThroat] = BsplineGeometry3(0, 'throat', nozzle.geometry.bSpline.knots, nozzle.geometry.bSpline.coefs);        
	    else
	        error('Only B-splines of degree 2 and 3 are supported');
	    end

	    if(xThroat ~= nozzle.geometry.xThroat)
	        fprintf('throat size/location changed with spline parameterization\n');
	    end
	    nozzle.geometry.xThroat = xThroat;
	    nozzle.geometry.xApparentThroat = nozzle.geometry.xThroat; % initialize apparent throat location
	    nozzle.throat.A = pi*yThroat^2;
        nozzle.inlet.A = pi*nozzle.geometry.bSpline.coefs(2,1)^2;
        nozzle.exit.A = pi*nozzle.geometry.bSpline.coefs(2,end)^2;
	    nozzle.geometry.Ainlet2Athroat = nozzle.inlet.A/nozzle.throat.A;
	    nozzle.geometry.Aexit2Athroat = nozzle.exit.A/nozzle.throat.A;    

	    % Make necessary functions for splined nozzle shape
	    A = @(x) BsplineGeometryMex(x,1,nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs');
	    dAdx = @(x) BsplineGeometryMex(x,2,nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs');
	    D = @(x) BsplineGeometryMex(x,3,nozzle.geometry.bSpline.knots,nozzle.geometry.bSpline.coefs');

	else % if nozzle shape is not a spline
	    A = @(x) nozzleGeometry(x,'A',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.shape);
	    dAdx = @(x) nozzleGeometry(x,'dAdx',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.shape);
	    D = @(x) nozzleGeometry(x,'D',nozzle.inlet.D,nozzle.geometry.length,nozzle.geometry.xThroat,nozzle.geometry.Ainlet2Athroat,nozzle.geometry.Aexit2Athroat,nozzle.geometry.shape);
	end

	
end