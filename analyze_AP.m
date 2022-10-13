function data = analyze_xs(filename,channels,channelnames,genotype,varargin)
%Finds the positions of nuclei and their concentration of dorsal.
%
%function data = analyze_xs(filename,channels,channelnames,genotype,varargin)
%
% This function first reads in an lsm-file of a zstack of embryo cross
% sections.  Next, it finds the periphery of the embryo and calculates the
% intensity of each (fluorsecent) color channel as you go around the
% periphery.  (This calculation is good for detection of mRNA expression
% patterns.)  Next, if a nuclear channel is present in the z-stack, this
% function locates the nuclei.  Further, if a nuclear-localized protein
% (such as dl or pMad) is present, the function then finds the nuclear
% concentration of that protein.  Finally, if an intronic probe is present,
% the function will find the nuclear dot intensity for each nucleus.  Note
% that neither the nuclear protein measurement nor the nuclear dot
% measurement will take place if no nuclear channel is present.
%
% Inputs:
% "filename": the name of the file, including the path to it
% "channels": 1-by-n vector describing what is in each channel, according
%	to the following codes:
%		1: nuclei
%		2: non-nuclear protein or mRNA
%		3: nuclear protein
%		4: intronic probe
%		5: none of these (and will not be processed)
%	So, for example, if in the first channel, you have nuclei, in the
%	second channel, you have a nuclear protein, in the third you have an
%	intronic probe, and in the fourth you have a DIC image, then your
%	vector should be: [1 3 4 5]
%
% Optional argument varargin can consist of these things, in this order:
%	* "yesplot": whether you want to plot the outcome.  Default, "false".  
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "ring_width": width of the annular ring.  That is, the distance in
%		microns into the embryo we take our data from. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "stage": what nuclear cycle is the embryo? it can be 10-14. Default,
%		14.
%		If this is not specified, but you still want to specify other 
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "nt": choice for number of bins in theta when detecting the embryo
%		periphery.  Default, 60. 
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "ns": choice for number of bins in "s" (pseudoarclength) when
%		measuring the intensity around the periphery.  Must be an even
%		number.  Default, 300.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
%
% Outputs:
% The function returns the structure 'data' with following the fields:
%
% filename: string vbl, name of the input file
% channels: the number code for what type of measurement is in each channel
% channeltypes: the type of each measurement in each channel
% channelnames: description of what's in each channel.
% genotype: string vbl, genotype of embryo
% metadata: a structure containing a bunch of the extra data that the user
%	doesn't often need to see or interact with (see below for explanation
%	of what goes in there).
% H: the number of y-pixels
% W: the number of x-pixels
% D: the number of z-slices
% Lbar: the average size of embryo (half-circumference) in microns
% s_mid: 1-by-n_channels vector of where the midline is predicted to be.
%	In most cases, the midline is the ventral midline, but in some (such as
%	pMad), it would be the dorsal midline.  Each channel gets an
%	opportunity to predict where the midline is.
% s: linspace(-1,1,301)', referring to the pseudo arclength.
% t: 301-by-n_channels array of the smoothened data as you go around the
%	periphery of the embryo.  Averaged in the z-direction.
% S: The pseudoarclength coordinates of each nucleus.  The z-direction is
%	ignored.
% R: The intensity of nuclear protein(s) for each nucleus.  Same length as
%	S, but has np columns, where np is the number of nuclear proteins.  The
%	identity of the nuclear protein in each column is in the same order as
%	specified in "channelnames".
% Std_R: the standard deviation of each measurement in "R".
% Intron: Similar to "R", but the intronic probe intensity for each
%	nucleus.
% A,B,M,mu,sig: gaussian fits for the nuclear proteins
% dA,dB,dM,dmu,dsig: length of 68% error-bars on above parameters
% gof: r-square goodness-of-fit value for the gaussian-fits.
% sV,sD,w: the ventral border, dorsal border, and widths of gene expression
%	patterns in the mRNA or intronic probes.  Each of these will be a
%	1-by-n_probe vector, where n_probe is the number of either mRNA or
%	intronic probe channels. The identity of the gene in each column is in
%	the same order as specified in "channelnames".
% dsV,dsD,dw: length of 68% error-bars on above parameters
% gof_gene: the r-square goodness-of-fit value for the gene expression fits
%
% The metadata field, which is a structure, has the following fields:
% lsminf1: series of metadata about the image taken; from 'lsmRead'
% lsminf2: more metadata; from 'lsminfo'
% scalings: 1x3 vector of the x,y,z scalings in microns per pixel
% rho: the nearest integer to the ratio of z scaling to xy scaling
% Yhatmax: the distance into the embryo (in pixels) the nuclear layer is
%	taken to be
% nt: number of points in theta where the embryo periphery is evaluated
% bg: the background levels of each of the channels.
% std_bg: the standard deviations of the background levels
% w: the average arclength, in pixels, of the embryo periphery
% arc: the arclength of each slice, in pixels
% L: the arclength of each slice in microns of the embryo's half-periphery
% Xp,Yp: cell array variables that contain the x and y coordinates of the
%	periphery of the embryo. Each element of these cell arrays corresponds
%	to a z-slice
% T,Raw: the intensity of the mRNA channel as you go around the
%	periphery of the embryo in quadrilaterals equally-spaced in periphery
%	pseudo arclength.  Each is a D-by-1 cell array.  Each element in the
%	cell arrays is a 301-by-n_channels array, where n_channels is the
%	number of non-"N/A" channels.
% X,Y,S: cell arrays containing the x,y and pseudoarclength coordinates of
%	each nucleus.  Each element of these cell arrays corresponds to a
%	z-slice. 
% Nuc,Std_nuc,Nuc_protein,Std_nuc_protein,Intron,Std_intron: cell array
%	variables that contain the intensity values and standard deviations of
%	nuclei, nuclear proteins, and intronic probe intensity, for each
%	nucleus.  Each element of these cell arrays corresponds to a z-slice.
% cint,cint68: the 95% and 68% confidence intervals on the gaussian fits to
%	the nuclear proteins: [A B M mu sig]
% B_gene: the background of gene expression fits.
% genes: the metadata for the gene fits
% introns: the metadata for the intronic probe fits





%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	yesplot = varargin{iArg}; else
 	yesplot = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	ring_width = varargin{iArg}; else
	ring_width = 18.36;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
 	stage = varargin{iArg}; else
 	stage = 14;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	nt = varargin{iArg}; else
 	nt = 60;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	ns = varargin{iArg}; else
 	ns = 500;
end%, iArg = iArg + 1;

if ~ischar(filename)
	error('wrong argument type')
end

%
% Reading in the lsm metadata and checking the "channels" input
%
lsminf2 = lsminfo(filename);
numch = lsminf2.NUMBER_OF_CHANNELS;
if length(channels) ~= numch
	fprintf('Number of channels in the z-stack: %d\n',numch)
	error('Make sure your "channels" input has the correct number of elements')
end

%
% Defining what is in each channel.
%
channeltype = cell(1,numch);
for i = 1:numch
	switch channels(i)
		case 1
			channeltype{i} = 'nuclei';
		case 2
			channeltype{i} = 'non-nuclear protein or mRNA';
		case 3
			channeltype{i} = 'nuclear protein';
		case 4
			channeltype{i} = 'intronic probe';
		case 5
			channeltype{i} = 'N/A';
		otherwise
			channels(i) = 5;
			channeltype{i} = 'N/A';
	end
end
nuc_ch = find(channels == 1);
if length(nuc_ch) > 1
	error('You can specify only one nuclear channel.')
end

%
% Reading in the lsm data
%
lsm =  lsmRead2(filename);
lsminf1 = lsm.lsm;
scalings = 1e6*[lsminf1.VoxelSizeX,lsminf1.VoxelSizeY,...
	lsminf1.VoxelSizeZ]; % microns/pixel
rho = round(scalings(3)/scalings(1)); % hope we don't have to round much.

IM = lsm.data;
clear lsm % to save memory.
IM(:,:,channels == 5,:) = [];
channels0 = channels;
channels(channels == 5) = [];
numch = size(IM,3);

%
% Extracting our region of interest.  There will be a rectangle
% that contains the embryo, outside of which all pixels will be
% zero.
%
if lsminf2.ScanInfo.USE_ROIS&&~lsminf2.ScanInfo.USE_REDUCED_MEMORY_ROIS
	I = sum(sum(IM,3),4);
	maxbw = logical(max(I)); j = find(maxbw); j1 = j(1); j2 = j(end);
	bwj = logical(I(:,j1)); i = find(bwj); i1 = i(1); i2 = i(end);
	
	IM = IM(i1:i2,j1:j2,:,:);
end

%
% Subtracting background.  We assume the background intensity level (i.e.,
% the true black level) for each channel is equal to the mode of
% intensities seen in a slice.
%
bg = zeros(numch,1);
bgsig = zeros(numch,1);
for i = 1:numch
	slice = double(IM(:,:,i,1));
	[n,x] = hist(slice(:),0:4095);
	[A,k] = max(n(1:end-1));
	bg(i) = x(k);
	X = x(max(k-4,1):k+4);
	Y = n(max(k-4,1):k+4);
	SIG = (X-bg(i))./sqrt(-2*log(Y/A));
	bgsig(i) = sqrt(meanDU(SIG'.^2));
	IM(:,:,i,:) = imsubtract(IM(:,:,i,:),bg(i));
end

D = lsminf1.DimensionZ; % depth
[H,W] = size(squeeze(IM(:,:,1,1))); % height, width
Yhatmax = round(ring_width/scalings(1)); % This is the width of the ring 
% (radius_outer - radius_inner) in pixels.  Default width: 18.36 microns.

%
% Middle remarks:
%
data.filename = filename;
data.channels = channels0;
data.channeltypes = channeltype;
data.channelnames = channelnames;
data.genotype = genotype;
data.metadata.lsminf1 = lsminf1;
data.metadata.lsminf2 = lsminf2;
data.metadata.scalings = scalings;
data.metadata.rho = rho;
data.metadata.Yhatmax = Yhatmax;
data.metadata.nt = nt;
data.metadata.bg = bg;
data.metadata.std_bg = bgsig;
data.H = H; 
data.W = W; 
data.D = D;


%
% First we find the periphery of the embryo and detect the intensity of
% each channel in quadrilaterals as you go around the periphery.  Next, we
% perform calculations on the nuclei/nuclear proteins/nuclear dots, if
% these things are present.  To do these things, we go through a loop that
% takes each slice individually.
%
h = 0.25;
arc = zeros(D,1);
T = cell(D,1); Raw = T;
Xp = T; Yp = Xp; 
		
X1 = cell(D,1); Y1 = X1; S = X1;
Nuc = X1; Std_nuc = X1;
Nuc_protein = X1; Std_nuc_protein = X1;
Intron = X1; Std_intron = X1;
for i = 1:D
	I = IM(:,:,:,i);
	
	%
	% border and perimiter.
	%
	[xp_in,yp_in] = borderFinder(I,h,yesplot,nt);
	[xp,yp] = borderFinder_AP(I,xp_in,yp_in,h,yesplot,nt);
% 	[xp,yp] = borderFinder_AP_old(I,h,yesplot,nt);
	Xp{i} = xp; Yp{i} = yp;
	arc(i) = sum(sqrt(diff(xp).^2 + diff(yp).^2));

	%
	% Peripheral intensities of each channel.
	%
	[t,raw,s] = domainMeas(I,xp,yp,Yhatmax,ns);
	T{i} = t;
	Raw{i} = raw;

	%
	% Now we detect the nuclei in each slice, given a nuclear channel.  
	%
	if ~isempty(nuc_ch)
		%
		% The function "find_nuclei" detects the nuclei (see the function
		% description for more info on how) and returns the pixel locations
		% of the nuclei, the locations of the centroids of the nuclei,
		% these locations in terms of pseudo-arclength, the perimeter of
		% the embryo, and a nuclear mask image.
		%
		[nucstats,xnuc,ynuc,snuc,w] = ...
			find_nuclei(I(:,:,nuc_ch),...
			xp,yp,scalings,Yhatmax,yesplot,stage);
		X1{i} = xnuc;
		Y1{i} = ynuc;
		S{i} = 2*snuc/w - 1;
		
		%
		% Obtaining the intensity of the nuclear protein channels, as well
		% as the nuclear channel itself, for each nucleus.
		%
		np_ch = find(channels == 3);
		[np,stdnp] = nuclearintensity(I(:,:,[nuc_ch np_ch]),nucstats);
		Nuc{i} = np(:,1);
		Std_nuc{i} = stdnp(:,1);
		if ~isempty(np_ch)
			Nuc_protein{i} = np(:,2:end);
			Std_nuc_protein{i} = stdnp(:,2:end);
		elseif i == D
			Nuc_protein = NaN; Std_nuc_protein = NaN;
		end
			
		
		%
		% Obtaining the intensity of any intronic probe channel.
		%
		v = channels == 4;
		if any(v)
			[Intron{i},Std_intron{i}] = intronicintensity(I(:,:,v),nucstats);
		else
			Intron = NaN; Std_intron = NaN;
		end

	elseif any(channels == 3)
		error('You must have a nuclear channel to have a nuclear protein channel!')
	elseif any(channels == 4)
		error('You must have a nuclear channel to have an intronic probe channel!')
	else
		X1 = NaN; Y1 = NaN; S = NaN;
		Nuc = NaN; Std_nuc = NaN;
		Nuc_protein = NaN; Std_nuc_protein = NaN;
		Intron = NaN; Std_intron = NaN;
	end
end

%
% Averaging the mRNA across all slices
%
t = T{1};
for i = 2:D
	t = t + T{i};
end
t = t/D;

%
% Transforming data on the nuclear protein into the scaled data
%
   if iscell(Nuc_protein)
  	p = cell2mat(Nuc_protein);
   	std_p = cell2mat(Std_nuc_protein);
   	nuc = cell2mat(Nuc);
   	r = p./repmat(nuc,1,size(p,2))*median(nuc);
   	std_nuc = cell2mat(Std_nuc);
   	stdev = r.*sqrt((std_p./p).^2 + ...
   		repmat((std_nuc./nuc).^2,1,size(p,2)));
   	
   	S1 = cell2mat(S);
   else
   	r = NaN;
   	stdev = NaN;
   	S1 = NaN;
   end

%
% Closing remarks
%
data.Lbar = mean(arc)/2*scalings(1);
data.s_mid = NaN(1,length(channels0));
data.s = s;
data.t = t;
data.S = S1;
data.R = r;
data.Std_R = stdev;

if iscell(Intron)
	data.Intron = cell2mat(Intron);
else
	data.Intron = Intron;
end


%
% Entering metadata.  This is good data for power-users to have around, but
% is not necessarily all that helpful under normal circumstances, so it
% remains largely hidden inside the substructure "metadata".
%
data.metadata.w = mean(arc);
data.metadata.arc = arc;
data.metadata.L = arc/2*scalings(1);
data.metadata.Xp = Xp;
data.metadata.Yp = Yp;
data.metadata.T = T;
data.metadata.Raw = Raw;
data.metadata.X = X1;
data.metadata.Y = Y1;
data.metadata.S = S;
data.metadata.Nuc = Nuc;
data.metadata.Std_nuc = Std_nuc;
data.metadata.Nuc_protein = Nuc_protein;
data.metadata.Std_nuc_protein = Std_nuc_protein;
data.metadata.Intron = Intron;
data.metadata.Std_intron = Std_intron;



%
% Below are placeholders
%
data.nucprotein_names = {};
data.A = NaN;
data.B = NaN;
data.M = NaN;
data.mu = NaN;
data.sig = NaN;
data.dA = NaN;
data.dB = NaN;
data.dM = NaN;
data.dmu = NaN;
data.dsig = NaN;
data.gof = NaN;

data.gene_names = {};
data.sV = NaN;
data.sD = NaN;
% data.sD2 = NaN;
% data.sD3 = NaN;
% data.sD4 = NaN;
% data.sD5 = NaN;
% data.sD6 = NaN;
% data.sD7 = NaN;
% data.sD8 = NaN;
% data.sD9 = NaN;
% data.sI = NaN;
% data.sF = NaN;
data.dsV = NaN;
data.dsD = NaN;
data.dsI = NaN;
data.dsF = NaN;
data.gof_gene = NaN;


save([data.filename(1:end-4),'_data.mat'],'data');



