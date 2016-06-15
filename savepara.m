function savepara(rivuletpara, outfilename)
	%SAVEPARA Summary of this function goes here
	%   Detailed explanation goes here

	if isempty(rivuletpara),
	    return;
	end
	outfilename(end-2:end) = [];
	outfilename = [outfilename 'txt']; 
	f = fopen(outfilename, 'wt');
	if f<=0,
	    error('Fail to open file to write');
	end	
	fprintf(f, 'plot: %0.2d percentage: %1.2f rewire: %1.1d gap: %5.0d ax_value: %1.0d dumpcheck: %1.0d\n, connectrate: %1.1f branchlen: %1.0d somagrowthcheck: %1.1f cleanercheck: %1.1f dtimageflag: %1.1f atmapflag: %1.1f\n, ignoreradiusflag: %1.1d prunetreeflag: %1.1d afmp: %1.2f speedatensorflag: %1.1d oofhmflag: %1.1f\n',...
		    	rivuletpara.plot_value, rivuletpara.percentage, rivuletpara.rewire, rivuletpara.gap, rivuletpara.ax_value,...
		    		rivuletpara.dumpcheck, rivuletpara.connectrate, rivuletpara.branchlen, rivuletpara.somagrowthcheck,... 
		    			rivuletpara.cleanercheck, rivuletpara.dtimageflag, rivuletpara.atmapflag, rivuletpara.ignoreradiusflag,...
		    				rivuletpara.prunetreeflag, rivuletpara.afmp, rivuletpara.speedastensorflag, rivuletpara.oofhmflag);
	fclose(f);

end