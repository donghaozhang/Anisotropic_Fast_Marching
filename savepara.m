function savepara(rivuletpara, outfilename)
	%SAVEPARA Summary of this function goes here
	%   Detailed explanation goes here

	if isempty(rivuletpara),
	    return;
	end

	f = fopen(outfilename, 'wt');
	if f<=0,
	    error('Fail to open file to write');
	end	
	fprintf(f, 'plot: %1.0d percentage: %0.2d rewire: %1.0d gap : %5.0d ax_value : %5.3f dumpcheck: %1.0d\n,...
		connectrate: %2.2d branchlen: %0.2d somagrowthcheck: %1.0d cleanercheck : %1.0d dtimageflag : %1.0d atmapflag: %1.0d\n,...
			connectrate: %2.2d branchlen: %0.2d somagrowthcheck: %1.0d cleanercheck : %5.0d dtimageflag : %5.3f atmapflag: %1.0d\n,...
		    	rivuletpara.plotvalue, rivuletpara.percentage, rivuletpara.rewire, rivuletpara.gap, rivuletpara.ax_value,...
		    		rivuletpara.plotvalue, rivuletpara.percentage, rivuletpara.rewire, rivuletpara.gap, rivuletpara.ax_value,...
		    			rivuletpara.dumpcheck, rivuletpara.connectr);
	fclose(f);

end