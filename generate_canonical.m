function [s,t,sD,sV,s_peak] = generate_canonical(soln,genename,gene_type,rbr,midline_channel)
% generates canonical profiles
%
%[s,t,sD,sV,s_peak] = generate_canonical(soln,genename,gene_type,rbr,midline_channel)
%
% This function takes 
%
% "soln": output of "analyze_xs.m" and probably "fit_gaussian.m" too.
% "genename": name of the gene (or protein) that you are attempting to make
%	the canonical profile for.
% "gene_type": There are three types of genes: (1) those whose expression
%	domains include the ventral midline, (2) those whose expression domains
%	include the dorsal midline, and (3) lateral genes.  The value of
%	"gene_type" is therefore either 'ventral', 'dorsal', or 'lateral', for
%	those three cases, respectively.
% "midline_channel": each channel has the opportunity to have an estimate
%	of the midline.  This index will tell us which channel to use as the
%	estimate for the ventral midline.

if ~exist('midline_channel','var') || isempty(midline_channel)
	midline_channel = NaN;
end

nLSM = length(soln);

ns = 300;
s = linspace(-1,1,ns+1)';
T = zeros(ns+1,nLSM);
% t1 = zeros(ns/2+1,nLSM);
% t2 = t1;
	
for i = 1:nLSM
	data = soln(i);
		
	%
	% Determining what to do with the midline.
	%
	if isnan(midline_channel) || isnan(data.s_mid(midline_channel))
		if all(isnan(data.s_mid))
			midline_channel = 1;
			data.s_mid = zeros(size(data.s_mid));
			error('Midline must be found first.  To find manually, use "plot_embryo".')
		else
			mc = find(~isnan(data.s_mid));
			midline_channel = mc(1); % default to first non-NaN channel
		end
	end
	
	
	s_mid = data.s_mid(midline_channel); % ventral midline
	sDM = s_mid + 1; sDM = mod(sDM+1,2) - 1; % dorsal midline
	[~,iDM] = min(abs(s-sDM));
	
	%
	% Determining what channel the gene (or protein) is in.
	%
	channelnames = data.channelnames;
	k = false(length(channelnames),1);
	for j = 1:length(channelnames)
		cname = channelnames{j};
		cname2 = str2cell(cname,',');  % split channelname in case mult
		% genes present.
		
		k(j) = ~isempty(strmatch(genename,cname2,'exact')); % this will be
		% logical true iff your genename appears exactly in the list of
		% genes in channel "j" (remember, each channel can potentially have
		% multiple genes/proteins).
	end
	if sum(k) ~= 1 % the gene must appear exactly once.
		%
		% if soln(i) does not contain your gene, or if for some reason your
		% gene appears twice in the same embryo, we stick NaN's in the ith 
		% column.
		%
		T(:,i) = NaN(ns+1,1);
	else
		channels = data.channels;
		channel = channels(k);
	
		switch channel
			case 1
				error('No canonical profiles for nuclei')
			case 2 % mRNA (or non-nuclear protein)
				%
				% When finding which column of "t" corresponds to our gene
				% or protein, we have to exclude any channels that are N/A
				% (usually brightfield channels)
				%
				v1 = false(size(channels));
				v1(k) = true;
				v1(channels == 5) = [];
				s1 = data.s;
				t1 = data.t(:,v1); 
				t = interp1(s1,t1,s);
				
				
				
			case 3 % nuclear protein
				%
				% There may be multiple nuclear proteins, so we have to
				% discover which-th one of them is the one of interest.
				%
				k_np = find(channels == 3);
				k_poi = find(k);
				whichnp = find(k_np == k_poi);
				
				
				%
				% Here we have to transform the nuclear protein
				% distribution, which is defined at discrete points along s
				% \in (-1,1], to be a smooth function defined at
				% equally-spaced points in "s".
				%
				S = data.S;
				R = data.R(:,whichnp);				
				[S,isort] = sort(S);
				R = R(isort);
				
				p = 50;
				S1 = [S(end-p+1:end)-2;S;S(1:p)+2]; % periodic extension
				R1 = [R(end-p+1:end);R;R(1:p)]; % periodic extension
				
				R2 = smooth(S1,R1,floor(p/2));
				[S2,R3] = repeat_remove(S1,R2); % necc to avoid duplicates
				t = interp1(S2,R3,s);
				
			case 4 % intronic probes
				%
				% There may be multiple nuclear proteins, so we have to
				% discover which-th one of them is the one of interest.
				%
				k_intr = find(channels == 4);
				k_goi = find(k);
				whichintr = find(k_np == k_goi);
				
				%
				% Intronic probes are even more complicated to get into a
				% smooth mesh, but we have already written a function for
				% that.
				%
				S = data.S;
				R = data.Intron(:,whichintr);
				t = smooth_intron(S,R,s);
			case 5
				error('No canonical profiles for N/A channel')
			otherwise
				error('No canonical profiles for undefined channel')
		end
		
		%
		% Aligning according to the midline.
		%
		t = t(1:end-1);
		t = [t(iDM:end);t(1:iDM)];
		T(:,i) = t;
		
% 		t1(:,i) = flipud(t(1:ns/2+1));
% 		t2(:,i) = t(ns/2+1:end);
	end
end

%
% Averaging all embryos together into one canonical peak
%
s = linspace(0,1,ns/2+1)';
t = meanDU(T,2);
t = subtrbkgrnd(t,rbr,true); % subtracting background
t1 = flipud(t(1:ns/2+1));
t2 = t(ns/2+1:end);
t = [t1 t2];
t = meanDU(t,2); % the "meanDU" function ignores NaNs.
if all(isnan(t))
	error('You may want to check if your gene/protein is in the input structure')
end
t = mDU(t); % normalizing between zero and one

%
% Plotting the peak so the user can interact and click on where the peak
% stops.  The canonical peak works best when the regions outside the peak
% are set to zero.
%
fignum = figure;
plot(s,t)
axis([0 1 0 1])
hold on
str1 = sprintf('(1) Left-click on the location(s) where the peak becomes zero.\n');
str1a = sprintf('\t\t(One on either side of the peak.)\n');
str2 = sprintf('(2) If you mess up, hit "d" to remove the last point.\n');
str3 = sprintf('(3) When you are done, hit "a".\n');

if strcmp(gene_type,'dorsal')
	title([str1 str2 str3])
	npts = 1;
	X = getpoint(npts);
	[i1,i1] = roundx(X(end),s);
	t(1:i1) = 0;
	
	s_offset = 1;
	s_peak = 1;
	sV = spline(t(i1+1:end),s(i1+1:end),0.5);
	sD = 1;
elseif strcmp(gene_type,'ventral')
	title([str1 str2 str3])	
	npts = 1;
	X = getpoint(npts);
	[i1,i1] = roundx(X(end),s);
	t(i1:end) = 0;
	
	s_offset = 0;
	s_peak = 0;
	sV = 0;
	sD = spline(t(1:i1-1),s(1:i1-1),0.5);
else % 'lateral', the default
	title([str1 str1a str2 str3])
	npts = 2;
	X = getpoint(npts);
	X = X(end-1:end); X = sort(X);
	[i1,i1] = roundx(X(1),s);
	[i2,i2] = roundx(X(2),s);
	t(1:i1) = 0;
	t(i2:end) = 0;
	
	[~,imax] = max(t);
	s_offset = s(imax);
	s_peak = s(imax);
	sV = spline(t(i1+1:imax),s(i1+1:imax),0.5);
	sD = spline(t(imax:i2-1),s(imax:i2-1),0.5);
end
close(fignum)
figure
t = mDU(t);
plot(s,t)
title(['Final canonical profile for ',genename])
save([genename,'avg'],'s','t','sD','sV','s_offset','s_peak','rbr')


% ------------- subfunction to get points off the graph ----------%
function X = getpoint(npts)

X = []; Y = []; h = [];
while true
	[x,y,button] = ginput(1);
	if button == ('a'+0) || button == 2 % we get out if points have already been taken
		if length(X) >= npts
			X = X(end-npts+1:end);
			break
		else
			str1 = sprintf('You need to pick %i point(s) before quitting.\n(Hit "enter" to continue)',npts);
			msgbox(str1,'Oops!','modal')
		end
	elseif button == ('d'+0) || button == 3 % we remove the last point
		if ~isempty(X)
			delete(h(end)); h(end) = [];
			X(end) = []; Y(end) = [];
		else
			msgbox('You can''t remove a point that''s not there!\n(Hit "enter" to continue)','Oops!','modal')
		end
	elseif button == 1 % we plot and store the point
		if y < 0, y = 0; end
		h1 = plot(x,y,'* r');
		h = [h;h1];
		X = [X;x]; Y = [Y;y];
% 		if exist('h1','var')
% 			delete(h1);
% 		end
	end
end

