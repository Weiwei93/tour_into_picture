function [corX, corY] = compute_outCorner(x1,y1,x2,y2,borderX,borderY)
% y = mx+b
m = (y1-y2)/(x1-x2);
b = y1-m*x1;

tempX = (borderY-b)/m;
tempY = m*borderX+b;

if(abs(tempX-x1) < abs(borderX-x1))
    corX = borderX;
    corY = tempY;
else
    corX = tempX;
    corY = borderY;
end