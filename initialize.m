function [J,Z,X]=initialize()
% This function is designed to independently assign motor positions and
% filament positions (randomly) within the hexagon.
% 01/26/14
% Callie J Miller
% USES THE FOLLOWING FUNCTIONS: ppolygon.m
% UPDATED 02/12/16: To implement potential force at boundary
global M N L x0 y0 xpol ypol Np hrad vrad distSide POLY maxX pointxM minX pointxm maxY pointyM minY pointym

%Pre-Allocate
J=zeros(2,M);
Z=zeros(5,N);
X=zeros(4,M);
%set boundaries of the patch of cortex
hrad=2; %radius of patch
vrad=2;
x0=0; %x component of center
y0=0; %y component of center
[xpol,ypol,POLY]=ellipse(hrad,vrad,x0,y0,Np);
[maxX,pointxM] = max(xpol);         %//Added the max points the filaments inside of the polygon can go
[minX,pointxm] = min(xpol);
[maxY,pointyM] = max(ypol);
[minY,pointym] = min(ypol);

distSide=sqrt((xpol(1)-xpol(2))^2+(ypol(1)-ypol(2))^2); % length of one side


%Initialize filament positions
for i=1:N
    Z(1,i) = (xpol(pointxm) + (xpol(pointxM)-xpol(pointxm))*rand); %//
    Z(2,i) = (ypol(pointyM) + (ypol(pointym)-ypol(pointyM))*rand); %//
    inPlus=inpolygon(Z(1,i),Z(2,i),xpol,ypol);
    while inPlus==0
        Z(1,i) = (xpol(pointxm) + (xpol(pointxM)-xpol(pointxm))*rand); %//
        Z(2,i) = (ypol(pointyM) + (ypol(pointym)-ypol(pointyM))*rand);
        inPlus=inpolygon(Z(1,i),Z(2,i),xpol,ypol);
    end
    Z(3,i)=rand()*2*pi;
    Z(4,i)=Z(1,i)-L*cos(Z(3,i));
    Z(5,i)=Z(2,i)-L*sin(Z(3,i));
    inMinus=inpolygon(Z(4,i),Z(5,i),xpol,ypol);
    while inMinus==0 % if minus ends aren't in hexagon, but the first while 
        % has already proven that the plus ends are inside the hexagon, 
        % define a new angle of orientation for the filaments, and then redefine a new minus end.
        Z(3,i) = rand()*2*pi;
        Z(4,i) = Z(1,i)-L*cos(Z(3,i));
        Z(5,i) = Z(2,i)-L*sin(Z(3,i));
        inMinus = inpolygon(Z(4,i),Z(5,i),xpol,ypol);
    end
end

%Initialize motor positions
for j = 1:M
  X(1,j) = (xpol(pointxm) + (xpol(pointxM)-xpol(pointxm))*rand);   %//
  X(2,j) = (ypol(pointyM) + (ypol(pointym)-ypol(pointyM))*rand);   %//
  inMotor = inpolygon(X(1,j),X(2,j),xpol,ypol);
  while inMotor==0
      X(1,j) = (xpol(pointxm) + (xpol(pointxM)-xpol(pointxm))*rand);   %//
      X(2,j) = (ypol(pointyM) + (ypol(pointym)-ypol(pointyM))*rand);   %//
      inMotor = inpolygon(X(1,j),X(2,j),xpol,ypol);
  end
end
X(3,:)=X(1,:); %Initially unattached
X(4,:)=X(2,:); %so same x and y position for L and R legs
