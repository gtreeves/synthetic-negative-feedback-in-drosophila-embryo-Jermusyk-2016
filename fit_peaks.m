function soln_out = fit_peaks(soln,varargin)
%Fits part of a color channel to a stereotypical peak of "gene".
%
%function soln_out = fitpeaks(soln,midline_channel,varargin)
%
% This function accepts as input the structure "soln", which is the output
% of the function "analyze_xs".  That structure contains intensity values
% as you go around the periphery of the embryo.  This function will fit
% peaks of those intensity values to stereotypical gene expression
% patterns.  For now, this only works for mature, nc 14 gene expression
% patterns.  The reason why is because there is no way to go from early,
% derepressed sog expresssion to completely repressed sog just by
% amplitude, offset, and stretch transformations.
%
% This function also fits intronic probe intensities, after performing some
% smoothing procedures on them.  The reason why we have to smooth them is
% because intronic probe intensities are inherently salt-and-pepper.
%
% This function can also fit profiles of nuclear proteins, after performing
% some smoothing procedures on them (and if you have a canonical profile
% for said nuclear protein).
%
% "soln": output of "analyze_xs.m".
%
% Optional argument varargin can consist of the following things:
%	* "yesplot": if logical "true", then a lot of plots will be
%		generated.  If false, plotting supressed.  Default, false.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "bkgrndthresh": fraction of maximum height where we say the
%		background begins.  Default, 0.15.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "peakthresh": fraction of maximum height of each individual peak
%		where we say the peak ends.  Default, 0.1.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% "soln_out": same as the input, "soln", but now the metadata under
%	soln_out.metadata.genes is filled in, as well as values for sV,sD, etc
%	under soln_out.


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	yesplot = varargin{iArg}; else
	yesplot = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	bkgrndthresh = varargin{iArg}; else
	bkgrndthresh = 0.15;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	peakthresh = varargin{iArg}; else
	peakthresh = 0.1;
end%, iArg = iArg + 1;

nstruct = length(soln);
soln_out = soln;

%
% Looping through each embryo in "soln"
%
for ii = 1:nstruct

	soln_out(ii).metadata.genes.bkgrndthresh = bkgrndthresh;
	soln_out(ii).metadata.genes.peakthresh = peakthresh;
	
	data = soln(ii);
	s_mid1 = data.s_mid;
	genotype = data.genotype;
	channels = data.channels;
	channels0 = channels;
	channels(channels0 == 5) = []; % removing the "N/A" channels
	s_mid1(channels0 == 5) = [];
	k_mrna = find(channels == 2); % "2" is the designation for mRNA
	k_np = find(channels == 3); % "3" is the designation for nuc protein
	k_intr = find(channels == 4); % "4" is the designation for intronic
	k_ch = sort([k_intr k_np k_mrna]);
	n_ch = length(k_ch); 
	
	%
	% Our data.
	%
	s = data.s;
	T = data.t;
	
	%
	% If we have intronic probes, we'll call a subfunction to smooth out
	% the salt-and-pepper nature of the intronic probe so we can fit those
	% profiles as well.  Once we do that, we replace our periphery
	% measurement of the intronic channels with the smoothed intron
	% version.
	%
	T_intr = data.Intron;
	if ~all(isnan(T_intr(:)))
		t_intr = smooth_intron(data.S,T_intr,s);
		T(:,k_intr) = t_intr;
	end
	
	%
	% If we have nuclear proteins, we'll call a subfunction to put the
	% profile on the same evenly-spaced mesh that mRNA is on.  Once we do
	% that, we replace our periphery measurement of the nuc protein
	% channels with the smoothed version.
	%
	S = data.S;
	R = data.R;
	if ~all(isnan(R))
		[S,isort] = sort(S);
		R = R(isort,:);
		p = 50;
		S1 = [S(end-p+1:end)-2;S;S(1:p)+2]; % periodic extension
		R1 = [R(end-p+1:end,:);R;R(1:p,:)]; % periodic extension
		
		R2 = zeros(size(R1));
		for k = 1:length(k_np)
			R2(:,k) = smooth(S1,R1(:,k),floor(p/2));
		end
		[S2,R3] = repeat_remove(S1,R2); % necc to avoid duplicates
		t_nucprot = interp1(S2,R3,s);
		try
		T(:,k_np) = t_nucprot;
		catch
			save
		end
	end
	
	%
	% This first loop over the different color channels gathers important
	% metadata about each channel.
	%
	params = cell(n_ch,1);
	S0 = params;
	Use_x0 = params;
	total_genes = 0;
	total_genenames = {};
	Rbr = zeros(n_ch,1);
	for j = 1:n_ch
		channelnames = data.channelnames;
		channelnames(channels0 == 5) = [];
		genenames = channelnames{k_ch(j)};
		
		%
		% Modifying our gene list to remove genes that have no expression
		% in the current genotype, and also noting if rho is one of the
		% genes.  We also look to see if any of our genes are repressed by
		% snail in a snail background.  Finally, we get an estimate of
		% where our peak should be.
		%
		genenames2 = str2cell(genenames,',');
		n_genes = length(genenames2);
		s0 = zeros(n_genes,1);
		use_x0 = true(n_genes,1);
		rbr0 = [];
		for k = 1:n_genes
			genename = genenames2{k};
			
			%
			% Getting an estimate of where our peak should be, based on
			% genename and genotype.
			%
			if strcmp(genename,'NaN')
				s0(k) = NaN;
			else
				load([genename,'avg'],'s_peak','rbr')
				rbr0 = [rbr0;rbr];
				s0(k) = s_peak;
			end
			if s0(k) == 1 || s0(k) == 0
				use_x0(k) = false;
			end
		end
		n_genes2 = n_genes;
		total_genes = total_genes + n_genes2;
		total_genenames = [total_genenames genenames2];
		S0{j} = s0;
		Use_x0{j} = use_x0;
		
		%
		% Building the parameters. The number and identity of the
		% parameters change depending on how many genes are in your color
		% channel and what type of genes they are (dorsal vs ventral vs
		% lateral).
		%
		C = cell(n_genes2,1);
		for k = 1:n_genes2
			kk = num2str(k);
			if use_x0(k)
				C{k} = ['A',kk,',delt',kk,',x0',kk,',']';
			else
				C{k} = ['A',kk,',delt',kk,',']';
			end
		end
		params{j} = char(C)';
		params{j}(end) = [];
		

		%
		% Subtracting background from raw curves.
		%
		Rbr(j) = max(rbr0);
		T(:,k_ch(j)) = subtrbkgrnd(T(:,k_ch(j)),Rbr(j));
	end
	
	%
	% Finding the ventral midline.  We had to loop through every channel to
	% gather enough data in order to reliably find the vm in each mRNA
	% channel. 
	%
	s_mid = find_midline(s,T(:,k_ch),S0,s_mid1(k_ch));
	kch5 = find(channels0 ~= 5); % to acct for removed "N/A" channels
	soln_out(ii).s_mid(kch5(k_ch)) = s_mid;
	
	%
	% Loop to fit each channel to canonical gene expression patterns.
	% Outputs are the ventral and dorsal borders as well as gene width.
	%
	sV = zeros(1,total_genes); sD = sV; w = sV;
	dsV = sV; dsD = sV; dw = sV; gof_gene = sV;
	count = 1;
	for j = 1:n_ch
		channelnames = data.channelnames;
		channelnames(channels0 == 5) = [];
		genenames = channelnames{k_ch(j)};
		genenames2 = str2cell(genenames,',');
		
		s0 = S0{j}; use_x0 = Use_x0{j}; n_genes2 = length(s0);
		t = T(:,k_ch(j));
		t = circShiftDU(t,s_mid(j)); % Realigning our data to the midline
				
		%
		% Splitting "t" in half
		%
		ns = length(s)-1;
		if isodd(ns)
			error('You must have an even number of points in x.')
		end
		npts = ns/2 + 1;
		s1 = linspace(0,1,npts)';
		t1 = t(npts:end); t2 = flipud(t(1:npts));
		
		%
		% Here we find the location and heights of our prospective peaks,
		% and store these in H1,I1 and H2,I2.
		%
		ds = 0.12;
		H1 = NaN(n_genes2,1); H2 = H1; I1 = H1; I2 = H1;
		for k = 1:n_genes2
			if s0(k) == 0
				H1(k) = t1(1);
				I1(k) = 1;
				H2(k) = t2(1);
				I2(k) = 1;
			elseif s0(k) == 1
				H1(k) = t1(end);
				I1(k) = npts;
				H2(k) = t2(end);
				I2(k) = npts;				
			else
				[k1,k1] = min(abs(s1-(s0(k)-ds)));
				[k2,k2] = min(abs(s1-(s0(k)+ds)));
				[H1(k),I1(k)] = max(t1(k1:k2));
				I1(k) = I1(k) + k1-1;
				[H2(k),I2(k)] = max(t2(k1:k2));
				I2(k) = I2(k) + k1-1;
			end
		end
		
		%
		% Calling a subfunction to do the fitting.  It's more like a
		% (3N+1)-parameter representation of our data, using the average,
		% stereotypical peaks of these genes as templates.  We do the fit
		% twice, once for each lateral half of the embryo.
		%
		n_params = sum(3*use_x0+2*~use_x0) + 1;
		k1 = isnan(H1);
		if ~all(k1)
			genenames3 = genenames2;
			genenames3(k1) = {'NaN'};
			genenames3 = cell2str(genenames3);
			fittypestr = ['genefit(''',genenames3,''',''',genotype,...
				''',x,B,',params{j},')'];
			f = fittype(fittypestr);
			[cfun1,cvals1,cint1,cint681,rsquare1] = fitelephant(t1,s1,H1,I1,f);
		else
			genenames3 = '';
			rsquare1 = NaN;
			cfun1 = NaN; 
			cvals1 = NaN(1,n_params); 
			cint1 = NaN(2,n_params); 
			cint681 = NaN(2,n_params); 
		end
		k2 = isnan(H2);
		if ~all(k2)
			genenames4 = genenames2;
			genenames4(k2) = {'NaN'};
			genenames4 = cell2str(genenames4);
			fittypestr = ['genefit(''',genenames4,''',''',genotype,...
				''',x,B,',params{j},')'];
			f = fittype(fittypestr);
			[cfun2,cvals2,cint2,cint682,rsquare2] = fitelephant(t2,s1,H2,I2,f);
		else
			genenames4 = '';
			rsquare2 = NaN; 
			cfun2 = NaN; 
			cvals2 = NaN(1,n_params);
			cint2 = NaN(2,n_params); 
			cint682 = NaN(2,n_params); 
		end
		
		%
		% Unpacking cvals1 and 2.
		%
		knan = strfind(genenames3,'NaN');
		A1 = cvals1(1:n_genes2); A1(knan) = NaN;
		B1 = cvals1(n_genes2+1);
		delt1 = cvals1(n_genes2+2:2*n_genes2+1); delt1(knan) = NaN;
		x01 = NaN(1,n_genes2);
		x01(use_x0) = cvals1(2*n_genes2+2:end); 
		x01(~use_x0) = s0(~use_x0);
		x01(knan) = NaN;
		
		knan = strfind(genenames4,'NaN');
		A2 = cvals2(1:n_genes2); A2(knan) = NaN;
		B2 = cvals2(n_genes2+1);
		delt2 = cvals2(n_genes2+2:2*n_genes2+1); delt2(knan) = NaN;
		x02 = NaN(1,n_genes2);
		x02(use_x0) = cvals2(2*n_genes2+2:end); 
		x02(~use_x0) = s0(~use_x0);
		x02(knan) = NaN;
		

		%
		% Plotting the result if asked for.
		%
% 		if yesplot
			plotfit(cfun1,cfun2,s1,t1,t2,ii,data.filename,genenames,genotype,rsquare1,rsquare2)
% 		end
		
		%
		% Calculating sV and sD, etc.
		%
		hh = 0.5;
		dcint1 = 0.5*diff(cint681);
		dcint2 = 0.5*diff(cint682);
		ddelt1 = dcint1(n_genes2+2:2*n_genes2+1);
		ddelt2 = dcint2(n_genes2+2:2*n_genes2+1);
		dx01(use_x0) = dcint1(2*n_genes2+2:end); 
		dx02(use_x0) = dcint1(2*n_genes2+2:end); 
		dx01(~use_x0) = 0;
		dx02(~use_x0) = 0;
		for k = 1:n_genes2
			
			%
			% Getting the borders of the canonical gene representation
			%
			genename = genenames2{k};
			[xV,xD,x_offset] = canonicalgeneborders(genename,genotype,hh);
			if isnan(xV), xV = 0; end
			if isnan(xD), xD = 1; end
			
			%
			% Calculating borders of this particlar gene expression peak
			%
			vloc1 = (xV-x_offset)*delt1(k) + x01(k);
			vloc2 = (xV-x_offset)*delt2(k) + x02(k);
			sV(count) = 0.5*(vloc1  + vloc2);
			
			dloc1 = (xD-x_offset)*delt1(k) + x01(k);
			dloc2 = (xD-x_offset)*delt2(k) + x02(k);
			sD(count) = 0.5*(dloc1 + dloc2);
			
			w(count) = sD(count) - sV(count);
			
			%
			% Calculating the errorbars on the parameters
			%
			dvloc1 = sqrt((xV-x_offset)^2*ddelt1(k).^2 + dx01(k).^2);
			dvloc2 = sqrt((xV-x_offset)^2*ddelt2(k).^2 + dx02(k).^2);
			dsV(count) = 0.5*sqrt(dvloc1^2 + dvloc2^2);
			
			ddloc1 = sqrt((xD-x_offset)^2*ddelt1(k).^2 + dx01(k).^2);
			ddloc2 = sqrt((xD-x_offset)^2*ddelt2(k).^2 + dx02(k).^2);
			dsD(count) = 0.5*sqrt(ddloc1^2 + ddloc2^2);
			
			dw = dsD + dsV;
			
			%
			% Goodness-of-fit (gof).  The gof will be the same for two
			% genes that are in the same color channel, but so we can keep
			% track of this for each gene, we'll record the values
			% separately anyway, even for two genes in the same color
			% channel.
			%
			gof_gene(count) = min(rsquare1,rsquare2);
			
			%
			% Advcancing the index by 1.
			%
			count = count + 1;
		end


		%
		% Metadata remarks
		%		
		soln_out(ii).metadata.genes.s0{j} = s0;
		soln_out(ii).metadata.genes.rbr(j) = Rbr(j);
		
		soln_out(ii).metadata.genes.A1{j} = A1;
		soln_out(ii).metadata.genes.B1{j} = B1;
		soln_out(ii).metadata.genes.x01{j} = x01;
		soln_out(ii).metadata.genes.delt1{j} = delt1;
		soln_out(ii).metadata.genes.gof1{j} = rsquare1;
		soln_out(ii).metadata.genes.cint1{j} = cint1;
		soln_out(ii).metadata.genes.cint681{j} = cint681;

		soln_out(ii).metadata.genes.A2{j} = A2;
		soln_out(ii).metadata.genes.B2{j} = B2;
		soln_out(ii).metadata.genes.x02{j} = x02;
		soln_out(ii).metadata.genes.delt2{j} = delt2;
		soln_out(ii).metadata.genes.gof2{j} = rsquare2;
		soln_out(ii).metadata.genes.cint2{j} = cint2;
		soln_out(ii).metadata.genes.cint682{j} = cint682;
	end
	
	%
	% Storing sV, sD, etc in our output
	%
	soln_out(ii).gene_names = total_genenames;
	soln_out(ii).sV = sV;
	soln_out(ii).sD = sD;
	soln_out(ii).w = w;
	soln_out(ii).dsV = dsV;
	soln_out(ii).dsD = dsD;
	soln_out(ii).dw = dw;
	soln_out(ii).gof_gene = gof_gene;
end

% ---------- subfunction to do the fitting -----------------
function [cfun,cvals,cint,cint68,r2] = fitelephant(t,s,A,I,f)

k = isnan(A);
A(k) = min(A(~k));
I(k) = 1;

%
% Upper and lower limits for the parameters (except x0)
%
A(A == 0) = 1e-2;
AL = 0.01*A; AU = 10*A;
BL = 0; BU = min(A); B = 0.5*BU;
delt = ones(1,length(A)); deltL = 0.5*delt; deltU = 2*delt;

%
% Determining what the values for x0 (initial guesses as well as upper and
% lower bounds) should be.  These change depending on whether we have
% lateral genes or not.  Lateral genes will have some intermediate value
% for x0 and will have upper and lower bounds.  Fully dorsal genes and
% fully ventral genes do not have an x0 for fitting purposes.
%
cnames = coeffnames(f);
nx = sum(strfindDU(cnames,'x'));
if nx == length(I)
	x0 = s(I);
else
	x0 = zeros(nx,1);
	count = 1;
	for i = 1:length(I)
		if I(i) ~=1 && I(i) ~= length(s) && ~k(i)
			x0(count) = s(I(i));
			count = count + 1;
		end
	end
end
x0L = x0 - 0.2; x0U = x0 + 0.2;

%
% Implementing our options
%
opts = fitoptions('Method','NonlinearLeastSquares',...
	'Startpoint',[A' B delt x0'],...
	'Lower',[AL' BL deltL x0L'],...
	'Upper',[AU' BU deltU x0U']);

%
% The actual fit, followed by unpacking the coefficient values, confidence
% intervals, and r2 value.
%
[cfun,gof] = fit(s,t,f,opts);
cvals = coeffvalues(cfun);
cint = confint(cfun);
cint68 = confint(cfun,(1-2*(1-normcdf(1,0,1))));
r2 = gof.rsquare;	
	





% ---------- subfunction to plot the fits -------------
function plotfit(cfun1,cfun2,s,t1,t2,i,filename,genenames,genotype,r21,r22)
figure(456)
set(gcf,'visib','off')
h = plot(s,[t1,cfun1(s),t2,cfun2(s)]);
set(h(1),'Color','c','Marker','o','Linestyle','none')
set(h(2),'Color','b','Linewidth',2)
set(h(3),'Color','m','Marker','o','Linestyle','none')
set(h(4),'Color','r','Linewidth',2)
title(['i = ',num2str(i),...
	',   genenames = ',genenames,...
	',   genotype = ',genotype,...
	',   r^{2}_1 = ',num2str(r21),...
	',   r^{2}_2 = ',num2str(r22)]);

v = strfind(genotype,'/');
genotype(v) = '_';
v = strfind(genotype,':');
genotype(v) = '-';

kslash = union(strfind(filename,'\'),strfind(filename,'/'));
filename(kslash) = '_';
filenameshort = filename(kslash(end-1)+1:end);

idx = num2strDU(i,3);
if ~(exist('Fittedpeaksimages','dir') == 7)
	mkdir Fittedpeaksimages
end
print(456,'-djpeg',['Fittedpeaksimages',filesep,filenameshort(1:end-4),'_',...
	genenames,'_',genotype,idx,'.jpg'])
close(456)

