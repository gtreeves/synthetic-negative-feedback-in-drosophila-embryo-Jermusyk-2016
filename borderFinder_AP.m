function [xp,yp] = borderFinder_AP(I,xp_in,yp_in,varargin)
%Finds the border that surrounds the cross-sectioned embryo.
%
%function [xp,yp] = borderFinder_AP(I,xp_in,yp_in,varargin)
%
% This function takes an image of am embryo cross-section and finds points
% on the outer edge of the embryo.  First, the image is gaussian
% filtered.  Then, taking 6-degree slices of the image (like a pizza),
% average intensity of the pixels within the slice is found as a function
% of distance from the center of the image.  The radius where the intensity
% drops to 25% maximal is where we say the edge of the embryo is.  We
% attempt some filtering of outliers at the end of the function. Single
% outliers can be detected by drastic changes in curvature.  Others can be
% corrected by fitting a circle to the embryo and removing those that are
% beyond two time the inner quartile range from the mean.
%
% "I": can be either a uint8 grayscale image, or string of filename.
%
% Optional argument varargin can consist of these things, in this order:
%	* "h": height of intensity cutoff.  Default, 0.25.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "yesplot": whether or not to plot the image and the boundary.
%		Default, false.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "nt": choice for number of bins for each quadrant of our elliptical
%		object.  Will determine length of outputs "xp" and "yp" (which will
%		be 4*nt+1).  Default, 15.
%
%
% "xp,yp": the points along the perimeter of the embryo.


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	h = varargin{iArg}; else
	h = 0.25;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	yesplot = varargin{iArg}; else
	yesplot = false;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	nt = varargin{iArg}; else
	nt = 60;
end%, iArg = iArg + 1;


dr = 0; % an arbitrary offset.

%
% Reading in "I"
%
if ischar(I)
	I = imread(I);
end
o = size(I,3);
if o > 1
	I = sum(double(I),3);
end

[m,n] = size(I);

%
% Gaussian filter, w/std 3.  Then, background subtract with a disk of
% radius 25 pixels.
%
sig = 3;
I = gaussFiltDU(I,sig);
I = imtophat(I,strel('disk',25));

%
% Fitting points to a ellipse.
%
[L,h1,xc,yc,phi] = ellipse_fit(xp_in(1:end-1),yp_in(1:end-1));
if phi > pi/2
	phi = phi - pi;
end

%
% foci of ellipse:
%
xf1 = xc + (sqrt(((.5*L)^2)-((.5*h1)^2)))*cos(phi);
yf1 = yc + (sqrt(((.5*L)^2)-((.5*h1)^2)))*sin(phi);
xf2 = xc - (sqrt(((.5*L)^2)-((.5*h1)^2)))*cos(phi);
yf2 = yc - (sqrt(((.5*L)^2)-((.5*h1)^2)))*sin(phi);

% =========================================================================
% Now we split up the embryo into four regions: near pole 1, between foci
% going from pole 1 to pole 2, near pole 2, and between foci going from
% pole 2 to pole 1.
% =========================================================================

%
% Region I: Near pole 1
%
x = (1:n) - xf1;  y = (1:m)' - yf1;
X = repmat(x,m,1); Y = repmat(y,1,n);

[Theta,R] = cart2pol(X,Y);
clear X Y

theta = linspace(-pi./2,pi./2,nt+1)';
theta = theta + phi;

r_out = borderFinder_quadrant(I,Theta,R,theta,h);

dtheta = theta(2) - theta(1);
xpI = r_out.*cos(theta(1:end-1)+dtheta/2) + xf1;
ypI = r_out.*sin(theta(1:end-1)+dtheta/2) + yf1;



%
% Region II: Between foci going from pole 1 to pole 2
%
x = (1:n) - xf1;  y = (1:m)' - yf1;
X = repmat(x,m,1); Y = repmat(y,1,n);

Xhat = -(X*cos(phi) + Y*sin(phi));
Yhat = -X*sin(phi) + Y*cos(phi);
% b = linspace(0,-sqrt((xf1-xf2)^2+(yf1-yf2)^2),nt+1)';
b = linspace(0,sqrt((xf1-xf2)^2+(yf1-yf2)^2),nt+1)';

r_out = borderFinder_quadrant(I,Xhat,Yhat,b,h);

db = b(2) - b(1);
xpII = -(b(1:end-1)+db/2)*cos(phi) - r_out*sin(phi) + xf1;
ypII = -(b(1:end-1)+db/2)*sin(phi) + r_out*cos(phi) + yf1;


%
% Region III: Near pole 2
%
x = (1:n) - xf2;  y = (1:m)' - yf2;
X = repmat(x,m,1); Y = repmat(y,1,n);

[Theta,R] = cart2pol(X,Y);
Theta(Theta<0) = Theta(Theta<0) + 2*pi;
clear X Y

theta = linspace(pi./2,3*pi./2,nt+1)';
theta = theta + phi;

r_out = borderFinder_quadrant(I,Theta,R,theta,h);

dtheta = theta(2) - theta(1);
xpIII = r_out.*cos(theta(1:end-1)+dtheta/2) + xf2;
ypIII = r_out.*sin(theta(1:end-1)+dtheta/2) + yf2;


%
% Region IV: Between foci going from pole 2 to pole 1
%
x = (1:n) - xf2;  y = (1:m)' - yf2;
X = repmat(x,m,1); Y = repmat(y,1,n);

Xhat = X*cos(phi) + Y*sin(phi);
Yhat = -(-X*sin(phi) + Y*cos(phi));
b = linspace(0,sqrt((xf1-xf2)^2+(yf1-yf2)^2),nt+1)';

r_out = borderFinder_quadrant(I,Xhat,Yhat,b,h);
r_out = -r_out;

db = b(2) - b(1);
xpIV = (b(1:end-1)+db/2)*cos(phi) - r_out*sin(phi) + xf2;
ypIV = (b(1:end-1)+db/2)*sin(phi) + r_out*cos(phi) + yf2;


% =========================================================================
% Now that we have the points in all four regions, we will combine them
% together and sort them into the right order.
% =========================================================================

xp = [xpI;xpII;xpIII;xpIV];
yp = [ypI;ypII;ypIII;ypIV];

theta = cart2pol(xp-xc,yp-yc);
[~,j] = sort(theta);
xp = xp(j);
yp = yp(j);

xp(end + 1) = xp(1);
yp(end + 1) = yp(1);


%
% Plotting, if you want.
%
if yesplot && ispc
	c1628(I)
	hold on
	plot(xp,yp,'r o -')
	plot(xf1,yf1,'*y')
	plot(xf2,yf2,'*y')
	plot(xc,yc,'*b')
	plot(xp_1,yp_1,'g')
	
end




% =========================================================================
% Subfunction to find the points on the periphery of the embryo in a given
% region
% =========================================================================

function r_out = borderFinder_quadrant(I,Theta,R,theta,h)

nt = length(theta) - 1;
r_out = zeros(nt,1);
for i = 1:nt
	
	%
	% I1 is a column vector of image intensity values that lie within our
	% current slice.  The col vec r1 coresponds to the radii of these
	% intensity values.  Then, we sort these vectors together such that the
	% values of r1 are ascending, and smooth.
	%
	v = Theta > theta(i) & Theta <= theta(i+1) & R >= 0;
	I1 = I(v);
	r0 = R(v);
	
	[r0,ir0] = sort(r0);
	r0 = round(r0);
	r = (r0(1):r0(end))'; nr = length(r);
	I1 = I1(ir0);
	npt = length(I1); % number of points for raw data.
	%n0 = floor(npt/nr); % properties of npt/nt.  This helps us determine the
	% shape and locations of "ones" in the matrix "M", which is the linear
	% transform from a string of individual pixel intensity values to those
	% averaged into bins.
	
	rowind = r0 - r0(1) + 1;
	colind = (1:npt)';
	M = sparse(rowind,colind,1,nr,npt); sumM = sum(M,2);
	Y0 = (M*double(I1))./sumM;
	Y0(isnan(Y0)) = 0;
	Y1 = imtophat(Y0,strel('line',100,90));
	
	%
	% Normalizing the intensities
	%
	y2 = mDU(Y1);
	i_out = find(y2(1:end-1) >= h & y2(2:end) < h);
	if ~isempty(i_out)
		i_out = i_out(1);
	else
		i_out = length(y2);
	end
	
	r_out(i) = r(i_out); % the radius of the edge of the embryo.

end




