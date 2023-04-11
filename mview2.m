function f = mview2(u,dyn,zoom,pauselength)
% mview2(u,dyn,zoom,pauselength)
%
% image sequence visuqlization
%
% dyn = intensity dynamics ([] means [min,max] like in imshow)
% z = zoom factor for display (either 1 number or [zx,zy])
% pauselength = Inf means that we start in pause mode
%
% type 'q' to stop
%      'space' to pause or play again
%      'r' to restart
%      left/right arrow to move frame by frame during pause
%
% author: Lionel Moisan                                   
% v1.0 (02/2022): first version from mview v1.1 (LM)

if nargin<2 || length(dyn)==0
    dyn = [min(u(:)),max(u(:))];
end
if nargin<3
    zoom = 1;
end
if nargin<4
    pauselength = 0;
end
if pauselength==Inf
    pauselength = 0;
    update = false;
else
    update = true;
end
last = true;
% to ensure correct key pressed handling
pauselength = max(0.05,pauselength);
if length(zoom)==1
    zoom = [zoom,zoom];
end
[ny,nx,nt] = size(u);
ny = ny*zoom(2);
nx = nx*zoom(1);
%screen = groot();
f = figure(...
    'Visible','off',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'DockControls','off',...
    'WindowStyle','normal',...
    'Colormap',gray(256));
im = image(kron(255*(u(:,:,1)-dyn(1))/(dyn(2)-dyn(1)),ones(zoom(2),zoom(1))),...
    'CDataMapping','direct',...
    'UserData',u,...
    'BusyAction','cancel',...
    'XData',0:nx-1,...
    'YData',0:ny-1);
axis 'image';
set(gca(),...
    'XLim',[0,nx-1],...
    'YLim',[0,ny-1],...
    'Visible','off',...
    'Units','pixels',...
    'ActivePositionProperty','position',...
    'Position',[1,1,nx,ny]);
set(f,...
    'KeypressFcn',@keypressed_callback,...
    'UserData',"");
%f.Position(1:2) = [min(100,screen.ScreenSize(3)-nx),min(screen.ScreenSize(4)-100-nx,screen.ScreenSize(4)-nx)];
f.Visible = 'on';
pause(0);
f.Position(3:4) = [nx,ny];
k = 1;
set(f,'Name',sprintf("image %d/%d",k,nt));
while true
    if pauselength==Inf
        pause()
    else
        pause(pauselength);
    end
    if length(f.UserData)>0
        key = f.UserData;
        if key=="q"
            close(f);
            break;
	elseif key=="space"
    	    update = ~update;
	elseif key=="r"
  	    k = 0;
	    last = false;
	elseif key=="leftarrow"
	    k = k-2;
	    last = false;
	elseif key=="rightarrow"
	    last = false;
	else
	    %fprintf("unrecognized key: %s\n",key);
	end
	f.UserData = '';
    end
    if update || ~last
        k = mod(k,nt)+1;
        im.CData = 255.*(kron(u(:,:,k),ones(zoom(2),zoom(1)))-dyn(1))/(dyn(2)-dyn(1));
        set(f,'Name',sprintf("image %d/%d",k,nt));
	last = true;
    end
end
end

function keypressed_callback(fig,event)
  fig.UserData = event.Key;
end
