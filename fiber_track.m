function [y,yraw,ok] = fiber_track(u,lambda,background)
% fiber tracking
%
% inputs:
%   u0: image sequence
%   lambda: curve smoothing parameter (e.g. 50)
%   background (optional): background average intensity + 1 std, estimated if not specified 
%
% outputs:
%   y: estimated fiber position 1st axis = column number, 2nd axis = image number
%   yraw (for control): same as y, but without the curve smoothing and outliers removal steps
%   ok: valid positions for yraw
%
% v1.0 (06/2022): original version (Lionel Moisan)
% v1.1 (06/2022): added peak robustness management (LM)

% parameters for peak robustness management
max_fiber_width = 10; % maximum fiber width, in pixels
fiber_relative_threshold = 0.5; % fiber defined by intensity > max_intensity * fiber_relative_threshold
max_smoothing_effect = 1; % if smoothing has more effect than this, class position as outlier

% filename = "20211112F9T2f-fixed4.tif";
% u = read_sequence(filename);
% lambda = 50;

% remove background (threshold is background level + 1 standard deviation of background level)
% could fail for very inhomogeneous backgrounds (in that case, specify threshold manually)
if nargin<3
    h = hist(u(:),1:ceil(max(u(:))));
    m = find(h==max(h));
    m = m(1);
    d = m-mean(u(u<=m));
    threshold = floor(m+d);
    fprintf("estimated threshold: %g\n",threshold);
    ut = max(0,u-threshold);
else %CB adjusts center of gravity for the fixed fiber
    threshold = zeros(size(u,1));
    for i=size(u,1)
        avg_ny = mean(u(i,:,:));
        threshold = u(i) - avg_ny(:);
    end
    ut = max(0,u-mean(threshold,1));
end
%ut = max(0,u-threshold);

% image smoothing (gaussian kernel parameters 31 and 17 may have to be adapted for different data)
[X,Y] = meshgrid(linspace(-3,3,31),linspace(-3,3,17));
ker = exp(-X.^2-Y.^2); ker = ker/sum(ker(:));
us = convn(ut,ker,'same'); 
[ny,nx,nt] = size(us);

% peak estimation (integer position)
[m0,y0] = max(us,[],1); 

% impose peak robustness
T = us>m0*fiber_relative_threshold;
[~,Y,~] = meshgrid(1:nx,1:ny,1:nt);
[~,y0m] = min(T.*Y+~T*1e10,[],1);
[~,y0p] = max(T.*Y-~T*1e10,[],1);
ok = squeeze(y0p-y0m<=max_fiber_width);


% sub-pixellic refinement of peak position using parabolic fitting
% note: we could add Fourier zoom here to increase precision
I = y0(:)'+ny*(0:nx*nt-1);
a = us(I-1);
b = us(I);
c = us(I+1);
dy = min(1,max(-1,0.5*(c-a)./(2*b-a-c)));
yraw = y0; yraw(:) = yraw(:)+dy';
yraw = squeeze(yraw); % remove 1st axis


% fiber smoothing and interpolation (use lambda parameter)
y = 0*yraw;
ker = exp(-linspace(-3,3,lambda).^2); ker = ker/sum(ker);
w = conv(yraw(:,1)*0+1,ker,'same'); % weight (for boundaries)
for t=1:nt
    s = yraw(:,t);
    for i=1:100
        s(ok(:,t)) = yraw(ok(:,t),t);
        s = conv(s,ker,'same')./w;
    end
    y(:,t) = s;
end

% remove new outliers and new interpolation
for t=1:nt
    ok(:,t) = ok(:,t) & abs(yraw(:,t)-y(:,t))<max_smoothing_effect;
    s = y(:,t);
    for i=1:100
        s(ok(:,t)) = yraw(ok(:,t),t);
        s = conv(s,ker,'same')./w;
    end
    y(:,t) = s;
end
   
fprintf("tracking complete.\n");


