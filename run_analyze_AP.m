function soln = run_analyze_AP(filename,channels,channelnames,genotype,outfilename)
%Runs "analyze_AP" on each lsm file in dir given in "pth".
%
%function soln = run_analyze_AP(filename,channels,channelnames,genotype,outfilename)
%
% "outfilename": Optional.  If you want the output save files to have a
%	name other than just "run_analyze_APCurrent.mat", it appends this
%	character variable onto the end of the name.

if ~exist('outfilename','var') || ~ischar(outfilename)
	outfilename = '';
end

if isdir(filename)
	pth = filename;
	filenames = readdir2(pth,'lsm');
else
	filenames = {filename};
end
nLSM = length(filenames);

%
% Looping through and extract data from each file.
%
j = 1; LE = {}; FE = {}; e = 1;
for i = 1:nLSM
	filename = filenames{i};
	
% 	try
		data = analyze_AP(filename,channels,channelnames,genotype);
		if any(channels == 3)
 			data = fit_gaussian(data);
		end
		if any(channels == 2) || any(channels == 4)
			data = fit_peaks(data,false);
		end
		save([filename(1:end-4),'_data.mat'],'data');
		soln(j,1) = data;
		save(['run_analyze_APCurrent',outfilename],'soln')
		disp(['j = ',num2str(j)])
		j = j + 1;
%  	catch
%  		lastE = lasterror;
%  		LE{end+1} = lastE;
%  		FE{end+1} = filename;
% 		fprintf('%s in %s\n',LE.message,filename)
%  		save(['run_analyze_APError',outfilename],'LE','FE')
%  		disp(['e = ',num2str(e)])
% 		e = e + 1;
% 	end

end








