function [ params ] = dakotaParseParams(FilNam)
	% Rick Feinrich & Victorien Menier, Feb 2016
	
	fid = fopen(FilNam);
	tline = fgets(fid);
	a = sscanf(tline, ' %d %s\n');
	num_vars = a(1);
	vars_text = char(a(2:end).');
	
	params.alt     = NaN;
	params.mach    = NaN;
	params.hInf    = NaN;
	params.Tstag   = NaN;
	params.Pstag   = NaN;
	params.t1      = NaN;
	params.t2      = NaN;
	params.t3      = NaN;
	params.rThroat = NaN;
	params.rExit   = NaN;
	params.kWall   = NaN;
	
	for ii = 1:num_vars
	    tline = fgets(fid);
	    a = sscanf(tline, ' %f %s\n');
	    var_i = a(1);
	    label_i = char(a(2:end).');

	    if(strcmp(label_i,'alt'))
	        params.alt  = var_i;
	    elseif(strcmp(label_i,'mach'))
	        params.mach = var_i;
	    elseif(strcmp(label_i,'hInf'))
	        params.hInf = var_i;
	    elseif(strcmp(label_i,'kWall'))
	        params.kWall = var_i;
	    elseif(strcmp(label_i,'TstagIn'))
	        params.Tstag = var_i;
	    elseif(strcmp(label_i,'PstagIn'))
	        params.Pstag = var_i;
	    elseif(strcmp(label_i,'t1'))
	        params.t1 = var_i;
	    elseif(strcmp(label_i,'t2'))
	        params.t2 = var_i;
	    elseif(strcmp(label_i,'t3'))
	        params.t3 = var_i;
	    elseif(strcmp(label_i,'rThroat'))
	        params.rThroat = var_i;
	    elseif(strcmp(label_i,'rExit'))
	        params.rExit = var_i;
	    else
	        disp('Unknown variable!');
	    end
	end
	fclose(fid);
	
	
end