% tested under MATLAB R2021a
% recommended to run script section by section, instead of all at once

% close all; clear all

% clear ycrop variable
if exist('ycrop')
    clear ycrop
end

% specify timelapse sequence here 
% define ycrop variable if you want to process only part of the lines of each frame
% example test files below, to recreate Figure 1E: 
filename = "20211112F8T140Hz live.tif";
%filename = "data_zebrafish/20211112F13f40Hz fixed.tif"; ycrop = 5:30;

% specify smoothing parameter here (e.g lambda = 50)
% (increase to obtain a smoother curve)
lambda = 50;

% vertical zoom used for display (change to 1 to avoid distortion)
zoom = 3; 

% set this variable to true to record video in file fibre.avi
record_video = false; 

% ------------- PROCESSING

% read timelapse
u0 = read_sequence(filename);

% crop some lines if necessary
if exist('ycrop')
    u0 = u0(ycrop,:,:);
end
    
% visualization
figure; imshow(mean(u0,3),[]);

%% fiber tracking
[y,yraw,ok] = fiber_track(u0,lambda);

%% controls
x = randi(size(u0,2));
figure; plot(squeeze(y(x,:))); pause(0.1); title(sprintf("fiber position across time for column %d",x));
pause(0.1);

t = randi(size(u0,3));
figure; plot(squeeze(y(:,t))); axis("equal"); pause(0.1); title(sprintf("fiber image %d",t));
pause(0.1);

h = floor(min(y(:))):0.001:ceil(max(y(:)));
figure; hist(y(:),h); pause(0.1); title("histogram of fiber positions (y)");

h = floor(min(yraw(ok(:)))):0.001:ceil(max(yraw(ok(:))));
figure; hist(yraw(ok(:)),h); pause(0.1); title("histogram of valid raw fiber positions (y)");

%% visualization with superimposition and zoom on y only
[ny,nx,nt] = size(u0);
clear F;
figure; 
for t=1:nt
    hold off;
    imshow(kron(u0(:,:,t),ones(zoom,1)),[]); hold on;
    %
    % visualization mode 1
    % plot((zoom+1)/2+(yraw(:,t)-1)*zoom,'g'); % to see the raw position before smoothing
    % plot(find(~ok(:,t)),(zoom+1)/2+(yraw(~ok(:,t),t)-1)*zoom,'b*'); % incorrect positions
    % plot((zoom+1)/2+(y(:,t)-1)*zoom,'r','linewidth',1); % estimated fiber
    %
    % visualization mode 2
    plot((zoom+1)/2+(y(:,t)-1)*zoom,'r','linewidth',1); % estimated fiber
    plot(find(~ok(:,t)),(zoom+1)/2+(y(~ok(:,t),t)-1)*zoom,'r+'); % interpolated positions
    title(sprintf("t=%d",t));
    %
    pause(0.1); % or pause() for keypressed wait
    %
    if record_video
        F(t) = getframe(gcf); pause(0.01);
    end
end
%% if requested, record display timelapse into file fibre.avi
if record_video
    writerObj = VideoWriter('fibre.avi');
    writerObj.FrameRate = 10;
    open(writerObj);
    for i=1:length(F)
        writeVideo(writerObj, F(i));
    end
    close(writerObj);
end
%% Figure 1. Plot RF position of live, paralyzed larva in rostrocaudal bins over time  

bin_size_live = 10;   
length_rc_axis = size(y,1);
multiple = length_rc_axis - mod(length_rc_axis, bin_size_live); %find multiple of rostrocaudal axis divided by bin size
timemultiple_length = multiple / bin_size_live;
reshape_multiple = reshape(y(1:multiple,:), bin_size_live, timemultiple_length, []); %reshape new binned double
bins_live = sum(reshape_multiple,1)/ bin_size_live; %take the mean over the 1st dimension
bins = reshape(bins_live, timemultiple_length, []);
   

%find which bin has the most movement (finding amplitude) to be able to
%plot this bin's position in the dorsoventral axis over time as an example
p = 0.194; %pixel size
micron_bins_live=bins * p; %convert to microns
sorted_bins_live=sort(micron_bins_live,2); %sort by rows to find the biggest amplitude
diff_bins_live=abs(bsxfun(@minus,sorted_bins_live, mean(sorted_bins_live))); %subtract to find the difference from the mean (absolute value)
max_bin_live=max(diff_bins_live(:,length(diff_bins_live)));

%search for max_bin_live value in diff_bins_live and pull out the row that it corresponds to it to plot
row_live=find(any(diff_bins_live == max_bin_live,2));
figure(1);
for i=row_live
    plot(micron_bins_live(i,:),'r'); hold on;
    title('RF dorsoventral position of live, paralyzed larva over time','FontSize', 20);
    xlabel('Number of frames','FontSize', 20)
    ylabel('Dorsoventral position (µm)','FontSize', 20)
    hold off; 
    set(gcf,'color','w');
    set(gca,'FontSize',20);
    axis tight
end
%% Figure 1. Repeat smoothing script for fixed larva and plot its dorsoventral position in bins over time

%use updated y variable from re-running the smoothing script, run this section INSTEAD of the section directly above in order to store the variables

%create bins for the example trace and find their average
bin_size_fixed = 10; 
length_rc_axis = size(y,1);
multiple = length_rc_axis - mod(length_rc_axis, bin_size_fixed); %find multiple of rostrocaudal axis divided by bin size
timemultiple_length = multiple / bin_size_fixed;
reshape_multiple = reshape(y(1:multiple,:), bin_size_fixed, timemultiple_length, []); %reshape new binned double
bins_fixed = sum(reshape_multiple,1)/ bin_size_fixed; %take the mean over the 1st dimension
bins = reshape(bins_fixed, timemultiple_length, []);
   
%find which bin has the least movement (finding amplitude) to be able to plot this bin's position in the dorsoventral axis over time
p = 0.194; %pixel size
micron_bins_fixed=bins * p; %convert to microns
sorted_bins_fixed=sort(micron_bins_fixed,2); %sort by rows
diff_bins_fixed=bsxfun(@minus,sorted_bins_fixed,sorted_bins_fixed(:,1));
min_bin_fixed=min(diff_bins_fixed(:,length(diff_bins_fixed)));

%search for max_bin_fixed value in diff_bins_fixed and pull out the row that it corresponds to it to plot
row_fixed=find(any(diff_bins_fixed == min_bin_fixed,2));
figure(2);
for i=row_fixed %the number of the row that had the most dorsoventral movement 
    plot(micron_bins_fixed(i,:),'k'); hold on;
    title('RF dorsoventral position of fixed larva over time','FontSize', 20);
    xlabel('Number of frames','FontSize', 20)
    ylabel('Dorsoventral position (µm)','FontSize', 20)
    hold off; 
    set(gcf,'color','w');
    set(gca,'FontSize',20);
    axis tight
end

%% Figure 1E1. Dorsoventral RF position of live, paralyzed larva vs fixed larva

figure(3);
for i=row_live 
    plot(micron_bins_live(i,:),'r','LineWidth',2); hold on;
    for i=row_fixed 
        micron_bins_adjusted = micron_bins_fixed * 1.24; %adjust view for figure
        plot(micron_bins_adjusted(i,:),'k','LineWidth',2); hold on;
    end
    set(gcf,'color','w');
    set(gca,'FontSize',20);
    legend('Live, paralyzed larva','Fixed larva', 'FontSize',20);
    %set(gca,'XColor','w','YColor','w','TickDir','out'); %comment out if you want axes to show
    ylabel(('Dorsoventral position'),'FontSize',20, 'Color', 'k')
    legend boxoff
    box off
    axis tight
end

%% Figure 1E2. Zoomed-in dorsoventral RF position of live, paralyzed larva vs fixed larva

figure(4);
for i=row_live 
    plot(micron_bins_live(i,(510:550)),'r-o','LineWidth',2); hold on; %pick any two timepoints along the graph that are 40 frames apart (imaging at 40Hz), showing a 1s zoomed-in trace
    for i=row_fixed
        micron_bins_adjusted = micron_bins_fixed * 1.28; %adjust view for figure
        plot(micron_bins_adjusted(i,(510:550)),'k-o','LineWidth',2); hold on;
    end
    set(gcf,'color','w');
    set(gca,'FontSize',20);
    legend('Live, paralyzed larva','Fixed larva', 'FontSize',20);
    %set(gca,'XColor','w','YColor','w','TickDir','out'); %comment out if you want axes to show
    ylabel(('Dorsoventral position'),'FontSize',20, 'Color', 'k')
    legend boxoff
    box off
    axis tight
end

%% Figure 1F1. Dorsoventral displacement probability density

%repeat the smoothing fiber function across multiple fish and store variables to reshape the matrix to have rostral, middle, caudal, and fixed columns

disp_bins = reshape(bins, [], 1);
micron_disp_bins = disp_bins * p;
mean_disp = mean(micron_disp_bins);
for j=1:size(micron_disp_bins)
    disp_from_mean = (micron_disp_bins(:) - mean_disp);
    mean_abs_dev = abs(disp_from_mean); 
    mean_abs_dev_std = sum(mean_abs_dev) / length(mean_abs_dev); %without removing outliers
    MAD = rmoutliers(mean_abs_dev, 'percentiles', [0 99]); %with removing outliers
    MADstd = sum(MAD) / length(MAD); 
end

% visualize data without removing outliers
figure;
histogram(mean_abs_dev);

% and with removing outliers
figure;
histogram(MAD);

%% Figure 1G. PCA 

filename = "20211112F8T140Hz live.tif";

% read images
clear u;
i = 0;
ok = true;
while ok
    i = i+1;
    try
        u(:,:,i) = double(imread(filename,i));
    catch
        i = i-1;
        ok = false;
    end
end
fprintf("%d images read from file\n%s\n",i,filename);

[ny,nx,nt] = size(u);

%% visualization
mview2(u,[],[1,5],0);

%% temporal smothing
sigma = 5; r = ceil(3*sigma); ker = exp(-(-r:r).^2/(2*sigma^2)); ker = ker/sum(ker);
um = convn(u,reshape(ker,1,1,[]),'same')./convn(ones(size(u)),reshape(ker,1,1,[]),'same');
mview2(um,[],[1,5],0);


%% visu mean/Fourier/std
figure; imshow(kron(mean(u,3),ones(5,1)),[]); pause(0.1); title('temporal mean')
figure; imshow(kron(fftshift(log(1+abs(fft2(mean(u,3))))),ones(5,1)),[]); pause(0.1); title('Fourier transform')
figure; imshow(kron([std(u,[],3);std(um,[],3)],ones(5,1)),[]); pause(0.1); title('temporal standard deviation (before/after temporal smoothing)');


%% spatial PCA (1 observation = 1 image) 
X = reshape(u,nx*ny,nt)';
X0 = mean(X,1);
[B,c] = pca(X-X0);
B = B';

[wcoeff,~,latent,~,explained]=pca(X-X0);

%% reconstruction with 20 components
n=1:20;
ur = reshape((X0+c(:,n)*B(n,:))',ny,nx,nt);
mview2([u;ur],[],[1,5],0);

%% same without the mean component
us = reshape((c(:,n)*B(n,:))',ny,nx,nt);
mview2(us,[],[1,5],0);

%% Figure 1G. visualization of the first 3 components
C = [];
for i=1:3
    d = reshape(B(i,:),ny,nx);
    C = [C;zeros(1,size(C,2));d/norm(d(:))];
end
C(ny+1:ny+1:end,:) = min(C(:));
figure; imshow(C,[-0.02, 0.02]); %adjust here for contrast for PCA in Fig 1


%% Movie 2. reconstruction along first 3 components
C = zeros(0,nx,nt);
for i = 1:3
    d = reshape((c(:,i)*B(i,:))',ny,nx,nt);
    if size(C,1)==0
        C = d;
    else
        C = [C;zeros(1,nx,nt);d];
    end
end
C(ny+1:ny+1:end,:,:) = min(C(:));
mview2(C,[],[1,2],0);


