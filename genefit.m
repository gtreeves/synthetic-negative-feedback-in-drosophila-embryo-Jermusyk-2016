function y = genefit(genenames,genotype,x,B,varargin)
%Transforms stereotypical peak of "gene".
%
%function y = genefit(genenames,genotype,x,B,varargin)
%
% This function loads preexisting data about the stereotypical shape and
% size of genes named "gene" (from files called "<gene>avg.mat") and
% transforms those data according to the amplitude, location, and
% "stretching factor" that the user provides ("A", "delt", and "x0", resp).
% In addition, the overall colorchannel will have a background, "B".
%
% Basically, the transformation is "z = (x - x0)/delt".
%
% There can be several different genes in a color channel, and this
% function will attempt to fit them all at the same time.  
%
% In the future, the following can be implemented.  First genes with
% multiple peaks, such as the gene "rho" has two separate domains, and as
% such, will be treated as two. Also could apply to eve, hb, etc. Second,
% sometimes "genotype" determines which data we will fit to.  For example,
% the genes sog, vnd, rho, and brk are all repressed by sna. Thus, in a sna
% background, these genes will extend to the ventral midline.  Neither of
% these have been implemented yet.
%
% "genenames": string, list of gene names (comma-separated) of genes you're
%	making an estimate of 
% "genotype": Not yet implemented, but possibly could load-in different
%	canonical profiles depending on the genotype.
% "x": the x-values where you'd like to estimate the peak
%
% "varargin" must be passed. The set of elements should be: 
%	"A": amplitude of the peak
%	"delt": stretching factor
%	"x0": location of the peak
% NOTE: only lateral genes will have an "x0" associated with them.
% 
% "y": your estimate.

%
% Input checks.
%
genenames2 = str2cell(genenames,',');
n_genes = length(genenames2);

y = zeros(size(x)) + B;
paramcount = 1;
for i = 1:n_genes
	genename = genenames2{i};
	if ~strcmp(genename,'NaN')

		%
		% Loadnig in data for your gene
		%
		matname = [genename,'avg'];
		load(matname,'s','t','s_offset')

		%
		% parameters
		%
		A = varargin{paramcount};
		delt = varargin{paramcount+1};
		if s_offset == 0 || s_offset == 1
			x0 = s_offset;
			paramcount = paramcount + 2;			
		else
			x0 = varargin{paramcount+2};
			paramcount = paramcount + 3;
		end
		
		%
		% The transformation and interpolation
		%
		z = (x - x0)./delt;
		y = y + A.*interp1(s-s_offset,t,z,'pchip',0);
	end
end









