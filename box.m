% get view box
close all;
I = imread('bbb.jpg');
%I = imresize(imread('bbb.jpg'), 0.5);   
imshow(I);
%[max_x, max_y, c] = size(I);
[max_y, max_x, c] = size(I);
[rx, ry] = ginput(2);

%% draw the rectangle
hold on;
irx = round([rx(1) rx(2) rx(2) rx(1) rx(1)]);
iry =  round([ry(1) ry(1) ry(2) ry(2) ry(1)]);
plot(irx,iry,'b'); 
hold off;

%% get vanishing point
while(1),
    [vx, vy, button] = ginput(1);
    
    if (isempty(button))
      break;
    end;
    vp_x = vx(1); vp_y=vy(1);

    % get the point on edge of image
%{
    cpx = zeros(4);
    cpy = [0 0 max_y max_y];
    cpx(1) = compute_corner(vp_x, vp_y, irx(1), iry(1), cpy(1)); 
    cpx(2) = compute_corner(vp_x, vp_y, irx(2), iry(2), cpy(2)); 
    cpx(3) = compute_corner(vp_x, vp_y, irx(3), iry(3), cpy(3)); 
    cpx(4) = compute_corner(vp_x, vp_y, irx(4), iry(4), cpy(4)); 
%}    
    cpx = [0 0 0 0];
    cpy = [0 0 0 0];
    [cpx(1), cpy(1)] = compute_corner(vp_x, vp_y, irx(1), iry(1), 1, 1); 
    [cpx(2), cpy(2)] = compute_corner(vp_x, vp_y, irx(2), iry(2), max_x, 1); 
    [cpx(3), cpy(3)] = compute_corner(vp_x, vp_y, irx(3), iry(3), max_x, max_y); 
    [cpx(4), cpy(4)] = compute_corner(vp_x, vp_y, irx(4), iry(4), 1, max_y); 
    
    % get the point on edge of image
    cpOutx = [0 0 0 0];
    cpOuty = [0 0 0 0];
    [cpOutx(1), cpOuty(1)] = compute_outCorner(vp_x, vp_y, irx(1), iry(1), 1, 1); 
    [cpOutx(2), cpOuty(2)] = compute_outCorner(vp_x, vp_y, irx(2), iry(2), max_x, 1); 
    [cpOutx(3), cpOuty(3)] = compute_outCorner(vp_x, vp_y, irx(3), iry(3), max_x, max_y); 
    [cpOutx(4), cpOuty(4)] = compute_outCorner(vp_x, vp_y, irx(4), iry(4), 1, max_y); 
    
    imshow(I);
    hold on;
    plot(irx,iry,'b'); 
    plot([vp_x cpx(1)],[vp_y cpy(1)],'r');
    plot([vp_x cpx(2)],[vp_y cpy(2)],'r');
    plot([vp_x cpx(3)],[vp_y cpy(3)],'r');
    plot([vp_x cpx(4)],[vp_y cpy(4)],'r');
    hold off;
end 


%% get depth of each face
focalLength = 400;

% left
if(cpOutx(1)==1)    
    i = cpOuty(1);
else
    i = cpy(1);
end

if(cpOutx(4)==1)
    j = cpOuty(4);
else
    j = cpy(4);
end

length1 = iry(4) - iry(1);
length2 = j - i;
ratio = length1/length2;
depthLeft = (focalLength-ratio*focalLength) / ratio; 

% right
if(cpOutx(2)==max_x)    
    i = cpOuty(2);
else
    i = cpy(2);
end

if(cpOutx(3)==max_x)
    j = cpOuty(3);
else
    j = cpy(3);
end

length1 = iry(3) - iry(2);
length2 = j - i;
ratio = length1/length2;
depthRight = (focalLength-ratio*focalLength) / ratio; 

% up
if(cpOuty(1)==1)    
    i = cpOutx(1);
else
    i = cpx(1);
end

if(cpOuty(2)==1)
    j = cpOutx(2);
else
    j = cpx(2);
end

length1 = irx(2) - irx(1);
length2 = j - i;
ratio = length1/length2;
depthUp = (focalLength-ratio*focalLength) / ratio; 

% bottom
if(cpOuty(4)==max_y)    
    i = cpOutx(4);
else
    i = cpx(4);
end

if(cpOuty(3)==max_y)
    j = cpOutx(3);
else
    j = cpx(3);
end

length1 = irx(3) - irx(4);
length2 = j - i;
ratio = length1/length2;
depthBottom = (focalLength-ratio*focalLength) / ratio;


%% crop each face of the box

% Middle face
mask = poly2mask([irx(1) irx(2) irx(3) irx(4)], [iry(1) iry(2) iry(3) iry(4)], size(I, 1), size(I, 2));
imMiddle = zeros(size(I,1), size(I,2), 3);
for i=1:size(I,1)
    for j=1:size(I,2)
        if(mask(i,j)==1)
            imMiddle(i,j,1) = double(I(i,j,1))/255;
            imMiddle(i,j,2) = double(I(i,j,2))/255;
            imMiddle(i,j,3) = double(I(i,j,3))/255;
        end
    end
end

figure(2), imshow(imMiddle);

% Left face
if(cpx(1)<0)
    if(cpx(4)>0)
        mask = poly2mask([cpx(1) irx(1) irx(4) cpx(4) 1], [cpy(1) iry(1) iry(4) cpy(4) size(I,1)], size(I, 1), size(I, 2));
    else
        mask = poly2mask([cpx(1) irx(1) irx(4) cpx(4)], [cpy(1) iry(1) iry(4) cpy(4)], size(I, 1), size(I, 2));
    end
else
    if(cpx(4)>0)
        mask = poly2mask([cpx(1) irx(1) irx(4) cpx(4) 1 1], [cpy(1) iry(1) iry(4) cpy(4) size(I,1) 1], size(I, 1), size(I, 2));
    else
        mask = poly2mask([cpx(1) irx(1) irx(4) cpx(4) 1 1], [cpy(1) iry(1) iry(4) cpy(4) 1 1], size(I, 1), size(I, 2));
    end
end
imLeft = zeros(size(mask,1), size(mask,2), 3);
for i=1:size(I,1)
    for j=1:size(I,2)
        if(mask(i,j)==1)
            imLeft(i,j,1) = double(I(i,j,1))/255;
            imLeft(i,j,2) = double(I(i,j,2))/255;
            imLeft(i,j,3) = double(I(i,j,3))/255;
        end
    end
end

% figure(3), imshow(imLeft);


% Right face
if(cpx(2)>size(I,2))
    if(cpx(3)>size(I,2))
        mask = poly2mask([cpx(2) cpx(3) irx(3) irx(2) size(I,2)], [cpy(2) cpy(3) iry(3) iry(2) 0], size(I, 1), size(I, 2));
    else
        mask = poly2mask([cpx(2) max_x cpx(3) irx(3) irx(2) size(I,2)], [cpy(2) max_y cpy(3) iry(3) iry(2) 0], size(I, 1), size(I, 2));
    end
else
    if(cpx(3)>size(I,2))
        mask = poly2mask([cpx(2) max_x cpx(3) irx(3) irx(2)], [cpy(2) 0 cpy(3) iry(3) iry(2)], size(I, 1), size(I, 2));
    else
        mask = poly2mask([cpx(2) max_x max_x cpx(3) irx(3) irx(2)], [cpy(2) 0 max_y cpy(3) iry(3) iry(2)], size(I, 1), size(I, 2));
    end
end
imRight = zeros(size(I,1), size(I,2), 3);
for i=1:size(I,1)
    for j=1:size(I,2)
        if(mask(i,j)==1)
            imRight(i,j,1) = double(I(i,j,1))/255;
            imRight(i,j,2) = double(I(i,j,2))/255;
            imRight(i,j,3) = double(I(i,j,3))/255;
        end
    end
end

% figure(4), imshow(imRight);


% Up face
if(cpx(2)<max_x)
    if(cpx(1)>1)
        mask = poly2mask([cpx(1) cpx(2) irx(2) irx(1)], [cpy(1) cpy(2) iry(2) iry(1)], size(I, 1), size(I, 2));
    else
        mask = poly2mask([cpx(1) 1 cpx(2) irx(2) irx(1)], [cpy(1) 1 cpy(2) iry(2) iry(1)], size(I, 1), size(I, 2));
    end
else
    if(cpx(1)>1)
        mask = poly2mask([cpx(1) max_x cpx(2) irx(2) irx(1)], [cpy(1) 1 cpy(2) iry(2) iry(1)], size(I, 1), size(I, 2));
    else
        mask = poly2mask([cpx(1) 1 max_x cpx(2) irx(2) irx(1)], [cpy(1) 1 1 cpy(2) iry(2) iry(1)], size(I, 1), size(I, 2));
    end
end
imUp = zeros(size(I,1), size(I,2), 3);
for i=1:size(I,1)
    for j=1:size(I,2)
        if(mask(i,j)==1)
            imUp(i,j,1) = double(I(i,j,1))/255;
            imUp(i,j,2) = double(I(i,j,2))/255;
            imUp(i,j,3) = double(I(i,j,3))/255;
        end
    end
end

% figure(5), imshow(imUp);


% Bottom face
if(cpx(3)<max_x)
    if(cpx(4)>1)
        mask = poly2mask([cpx(3) cpx(4) irx(4) irx(3)], [cpy(3) cpy(4) iry(4) iry(3)], size(I, 1), size(I, 2));
    else
        mask = poly2mask([cpx(3) 1 cpx(4) irx(4) irx(3)], [cpy(3) max_y cpy(4) iry(4) iry(3)], size(I, 1), size(I, 2));
    end
else
    if(cpx(4)>1)
        mask = poly2mask([cpx(3) max_x cpx(4) irx(4) irx(3)], [cpy(3) max_y cpy(4) iry(4) iry(3)], size(I, 1), size(I, 2));
    else
        mask = poly2mask([cpx(3) max_x 1 cpx(4) irx(4) irx(3)], [cpy(3) max_y max_y cpy(4) iry(4) iry(3)], size(I, 1), size(I, 2));
    end
end

imBottom = zeros(size(I,1), size(I,2), 3);
for i=1:size(I,1)
    for j=1:size(I,2)
        if(mask(i,j)==1)
            imBottom(i,j,1) = double(I(i,j,1))/255;
            imBottom(i,j,2) = double(I(i,j,2))/255;
            imBottom(i,j,3) = double(I(i,j,3))/255;
        end
    end
end

% figure(6), imshow(imBottom);



%% get texture maps
% homographies for each face

% Up face
pts1 = [cpx(1) cpx(2) irx(2) irx(1) ; cpy(1) cpy(2) iry(2) iry(1) ; 1 1 1 1];
%pts2 = [cpx(1) cpx(2) cpx(2) cpx(1) ; cpy(1) cpy(2) iry(2) iry(1) ; 1 1 1 1];
pts2 = [irx(1) irx(2) irx(2) irx(1) ; 1 1 depthUp depthUp ; 1 1 1 1];
H_Up = computeHomography(pts1, pts2);

T = maketform('projective', H_Up.');
imUpBlend = imtransform(imUp, T, 'XData',[irx(1) irx(2)],'YData',[0 depthUp]);

figure(3), imshow(imUpBlend);


% Bottom face
pts1 = [cpx(3) cpx(4) irx(4) irx(3) ; cpy(3) cpy(4) iry(4) iry(3) ; 1 1 1 1];
%pts2 = [cpx(3) cpx(4) cpx(4) cpx(3) ; cpy(3) cpy(4) iry(4) iry(3) ; 1 1 1 1];
pts2 = [irx(3) irx(4) irx(4) irx(3) ; depthBottom depthBottom 1 1 ; 1 1 1 1];
H_Bottom = computeHomography(pts1, pts2);

T = maketform('projective', H_Bottom.');
imBottomBlend = imtransform(imBottom, T, 'XData',[irx(4) irx(3)],'YData',[0 depthBottom]);

figure(4), imshow(imBottomBlend);


% Left face
pts1 = [cpOutx(1) irx(1) irx(4) cpOutx(1) ; cpOuty(1) iry(1) iry(4) cpOuty(4) ; 1 1 1 1];
%pts2 = [1 irx(1) irx(4) 1 ; 1 1 max_y max_y ; 1 1 1 1];
pts2 = [1 depthLeft depthLeft 1 ; iry(1) iry(1) iry(4) iry(4) ; 1 1 1 1];
H_Left = computeHomography(pts1, pts2);

T = maketform('projective', H_Left.');
imLeftBlend = imtransform(imLeft, T, 'XData',[0 depthLeft],'YData',[iry(1) iry(4)]);
figure(5), imshow(imLeftBlend);


% Right face
pts1 = [cpOutx(3) irx(3) irx(2) cpOutx(2) ; cpOuty(3) iry(3) iry(2) cpOuty(2) ; 1 1 1 1];
%pts2 = [max_x irx(3) irx(2) max_x ; max_y max_y 1 1 ; 1 1 1 1];
pts2 = [depthRight 1 1 depthRight ; iry(3) iry(3) iry(2) iry(2) ; 1 1 1 1];
H_Right = computeHomography(pts1, pts2);

T = maketform('projective', H_Right.');
imRightBlend = imtransform(imRight, T, 'XData',[0 depthRight],'YData',[iry(2) iry(3)]);
figure(6), imshow(imRightBlend);


% Middle face
imMiddleBlend = imcrop(imMiddle, [irx(1) iry(1) (irx(2)-irx(1)) (iry(4)-iry(1))]);
figure(7), imshow(imMiddleBlend);

%% define 3D surfaces
xmin = 0; xmax = irx(2)-irx(1);
ymin=0; ymax = iry(4)-iry(1);

up_x = [depthUp depthUp depthUp; 0 0 0];
up_y = [xmax (xmax+xmin)/2 xmin; xmax (xmin+xmax)/2 xmin];
up_z = [ymax ymax ymax; ymax ymax ymax];

bottom_x = [depthBottom depthBottom depthBottom; 0 0 0];
bottom_y = [xmax (xmax+xmin)/2 xmin; xmax (xmin+xmax)/2 xmin];
bottom_z = [ymin ymin ymin; ymin ymin ymin];

middle_x = [0 0 0; 0 0 0];
middle_y = [xmax (xmax+xmin)/2 xmin; xmax (xmin+xmax)/2 xmin];
middle_z = [ymax ymax ymax; ymin ymin ymin];

left_x = [0 depthLeft/2 depthLeft; 0 depthLeft/2 depthLeft];
left_y = [xmin xmin xmin; xmin xmin xmin];
left_z = [ymax ymax ymax; ymin ymin ymin];

right_x = [depthRight depthRight/2 0; depthRight depthRight/2 0];
right_y = [xmax xmax xmax; xmax xmax xmax];
right_z = [ymax ymax ymax; ymin ymin ymin];

% create the surface and texturemap it with a given image

view = figure('name','3DViewer: Directions[W-S-A-D] Zoom[Q-E] Exit[ESC]');
set(view,'windowkeypressfcn','set(gcbf,''Userdata'',get(gcbf,''CurrentCharacter''))') ;
set(view,'windowkeyreleasefcn','set(gcbf,''Userdata'','''')') ;
set(view,'Color','black')
hold on

warp(up_x, up_y, up_z, imUpBlend);
warp(bottom_x, bottom_y, bottom_z, imBottomBlend);
%warp(right_x, right_y, right_z, imLeftBlend );
warp(left_x, left_y, left_z,imRightBlend );
warp(middle_x, middle_y, middle_z, imMiddleBlend);

%% camera
% Some 3D magic...
axis equal;  % make X,Y,Z dimentions be equal
axis vis3d;  % freeze the scale for better rotations
axis off;    % turn off the stupid tick marks
camproj('perspective');  % make it a perspective projection

% set camera position
camx = 60;
camy = 145;
camz = 38.8;

% set camera target
tarx = 28.5;
tary = 112;
tarz = 117.5;

% set camera step
stepx = 0.05;
stepy = 0.05;
stepz = 0.05;

% set camera on ground
camup([0,0,1]);
campos([camx camy camz]);

key = 0;
while (~key),
    waitforbuttonpress;
    key = get(view, 'currentch');
    
    switch key
        case 'd'
            camdolly(-stepx,0,0,'fixtarget');
        case 'a'
            camdolly(stepx,0,0,'fixtarget');
        case 's'
            camdolly(0,stepy,0,'fixtarget');
        case 'w'
            camdolly(0,-stepy,0,'fixtarget');
        case 'q'
            camdolly(0,0,stepz,'fixtarget');
        case 'e'
            camdolly(0,0,-stepz,'fixtarget');
        case 28
            campan(-0.1,0);
        case 29
            campan(0.1,0);
        case 30
            campan(0,0.1);
        case 31
            campan(0,-0.1);
        case 13
            break;
        case 27
            break;
        case 'p'
            break;
    end
    
    key = 0;

    pause(.001);

    %campos([camx camy camz]);
    %camtarget([tarx tary tarz]);
    pos = campos;
    target = camtarget;
    
end;


