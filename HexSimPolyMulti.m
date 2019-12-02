function [X,J,Z,S]=HexSimPolyMulti(X,J,Z,Stiffness)
% This code (saved in /Users/CallieMiller/Desktop/Demetrius'
% Code/Simulation Codes/HexSimPoly.m) is used to quickly record the 3 data
% values for filaments (plus end location Z(1,i),Z(2,i) and the angle of
% the filament Z(3,i)) to later be read into ImageJ for quick plotting and
% analyzing. It is a similar version to the HexSim.m code in Demetrius'
% Code folder and Image Processing folder. It also includes a noise concept
% of filament polymerization.
% 01/26/14
% MODIFIED: 2/13/14 to include a threshold of stretching for the motor
% which is the same as the radius of searching for attachments, r.
% MODIFIED: 2/17/14 to include a .txt recording of ForceMag(i)- the sum of
% all motor induced forces on filament i
% Callie J Miller
% USES THE FOLLOWING FUNCTIONS: findfil2D.m, newMotorPos.m, ForceVertex.m
% (betaangle.m)
%
% MODIFIED: 2/26/14 to bias motor diffusion to one side of the hexagon
% MODIFIED: 3/6/14 to correct how I interpret motor diffusion
% MODIFIED: 3/25/14 to correct the time scaling for motor attach/detach.
% Was p0>rand*h, for example, corrected to p0*h>rand.
% UPDATED: 05/01/15 to correct for force calculations on filaments. Changed
% filament movement section of the code
% UPDATED: 02/12/16 to implement potential force at boundary

global runSims newFolder th BoundaryRadius Time M N L p0 p1 p2 h r v k k0 hrad vrad x0 y0 xpol ypol Zold Np t DeltaX distSide maxX pointxM minX pointxm maxY pointyM minY pointym
% Z- matrix giving the positional information for each filament
%   Z(1,i)=x-position of plus end
%   Z(2,i)=y-position of plus end
%   Z(3,i)=angle of filament
%   Z(4,i)=x-position of minus end
%   Z(5,i)=y-position of minus end
% J- number of filament motor is attached to
%   J(1,i)= L leg of motor: 0 if unattached
%   J(2,i)= R leg of motor: 0 if unattached
% X- matrix giving positional information for each motor
%   X(1,j)= x-position of L leg
%   X(2,j)= y-position of L leg
%   X(3,j)= x-position of R leg
%   X(4,j)= y-position of R leg
% N- number of filaments
% M- number of motors
% Time- time to run the simulation
d=0.02; % rate of diffusion
mu=0.001; % step size of diffusion
Dif=d/(d+p1); %Diffusion value for compare.

% frictional coefficients based on Hunt, Gittes, Howard paper
eta=1; %x10^(-15) Ns/um^2
di=0.008; %um, diameter of actin
p=L/di;
gperp=0.84;
gparl=-0.2;
grot=-0.662;
l_pr=(2*pi*eta*L)/(log(p)+gparl);
l_pp=(4*pi*eta*L)/(log(p)+gperp);
l_r=((1/3)*pi*eta*L^3)/(log(p)+grot);
perim=2*pi*hrad; % perimeter will need to change if it's an ellipse

% Defining vertex parameters
E=zeros(1,Np+1);
for n=1:Np
    if (n<15 || n>89 || (n>39 && n<65)) % Sides of 100 sided polygon
        E(n)=Stiffness; 
    else
        E(n)=100; 
    end
end


% Pre-allocating
S=zeros(N,M);
mxA=zeros(1,M);
myA=zeros(1,M);
mxB=zeros(1,M);
myB=zeros(1,M);
x1=zeros(1,M);
y1=zeros(1,M);
F=zeros(2,M);
FR=zeros(2,M);
ang=zeros(1,M);
count=zeros(1,N);
FB=zeros(2,N);
angFB=zeros(1,N);

for t = 1:Time
    %MOTOR ATTACHING
    for j=1:M
        if J(1,j)==0 && J(2,j)==0 % R and L legs of motor are unattached
            [count]=findfil2Dv2(Z,X(1,j),X(2,j));
            if sum(count)==0 || rand()<Dif || p1*h<rand() %if these things are true, motor will walk
                phi=rand()*2*pi;
                step=normrnd(mu,1)*h;
                X(1,j)=X(1,j)+step*cos(phi);
                X(2,j)=X(2,j)+step*sin(phi);
                
                inMotor = inpolygon(X(1,j),X(2,j),xpol,ypol);
                while inMotor==0
                    X(1,j) = (xpol(pointxm) + (xpol(pointxM)-xpol(pointxm))*rand);   %//
                    X(2,j) = (ypol(pointyM) + (ypol(pointym)-ypol(pointyM))*rand);   %//
                    inMotor = inpolygon(X(1,j),X(2,j),xpol,ypol);
                end
                X(3,j)=X(1,j);
                X(4,j)=X(2,j);
                J(1,j)=0;
                J(2,j)=0;
            elseif sum(count)==1 % then both legs of motor attach to 1 filament
                a=find(count);
                J(1,j)=a;
                J(2,j)=a;
                [X(1,j),X(2,j)]=newMotorPos(Z(1,a),Z(2,a),Z(4,a),Z(5,a),X(1,j),X(2,j));
                X(3,j)=X(1,j);
                X(4,j)=X(2,j);
            else
                a=find(count); %will find the index of the nonzero elements in count
                indicies=ceil(rand()*max(size(a))); % randomly pick one of these for L leg to bind
                fil=a(indicies);
                J(1,j)=fil;
                % Define new position for L motor once it attaches
                [X(1,j),X(2,j)]=newMotorPos(Z(1,fil),Z(2,fil),Z(4,fil),Z(5,fil),X(1,j),X(2,j));

                indicies2=ceil(rand()*max(size(a))); % randomly pick one of these for R leg to bind
                fil2=a(indicies2);
                while fil2==fil
                    indicies2=ceil(rand()*max(size(a))); % randomly pick one of these for R leg to bind
                    fil2=a(indicies2);
                end
                J(2,j)=fil2;
                % Define new position for R motor once it attaches
                [X(3,j),X(4,j)]=newMotorPos(Z(1,fil2),Z(2,fil2),Z(4,fil2),Z(5,fil2),X(3,j),X(4,j));
            end
        elseif J(1,j)==0 && J(2,j) ~=0 % L leg motor is unattached but R is attached
            [count]=findfil2Dv2(Z,X(1,j),X(2,j));
            if sum(count)==0 || p1*h<rand() % No filaments in the area or rate of attachment doesn't occur, then L leg of motor can't attach to any filament
                J(1,j)=0; % L motor stays unattached
                J(2,j)=J(2,j); % R motor stays attached to original filament
            else % L motor can attach to a filament
                a=find(count); %will find the index of the nonzero elements in count
                indicies=ceil(rand()*max(size(a))); % randomly pick one of these for L leg to bind
                fil=a(indicies);
                J(1,j)=fil;
                % Define new position for L motor once it attaches
                [X(1,j),X(2,j)]=newMotorPos(Z(1,fil),Z(2,fil),Z(4,fil),Z(5,fil),X(1,j),X(2,j));
            end
            if p0*h >rand() % R leg unattaches
                J(2,j)=0;
                X(3,j)=X(1,j);
                X(4,j)=X(2,j);
            end
        elseif J(1,j)~=0 && J(2,j)==0 % R leg motor is unattached but L is attached
            [count]=findfil2Dv2(Z,X(3,j),X(4,j));
            if sum(count)==0 || p1*h<rand() % then R leg of motor can't attach to any filament
                J(1,j)=J(1,j); % L motor stays attached to original filament
                J(2,j)=0; % R motor stays unattached
            else % R motor can attach to a filament
                a=find(count); %will find the index of the nonzero elements in count
                indicies=ceil(rand()*max(size(a))); % randomly pick one of these for R leg to bind
                fil=a(indicies);
                J(2,j)=fil;
                % Define new position for R motor once it attaches
                [X(3,j),X(4,j)]=newMotorPos(Z(1,fil),Z(2,fil),Z(4,fil),Z(5,fil),X(3,j),X(4,j));
            end
            if p0*h >rand() % L leg unattaches
                J(1,j)=0;
                X(1,j)=X(3,j);
                X(2,j)=X(4,j);
            end
        elseif J(1,j)==J(2,j) && J(1,j)~=0 % R and L legs on the same filament
            if rand() >0.5 % then R leg motor can attach to a second filament
                [count]=findfil2Dv2(Z,X(3,j),X(4,j));
                if p1*h > rand() && sum(count)~=0% R leg can attach if there's a filament in the area
                    a=find(count);
                    indicies=ceil(rand()*max(size(a)));
                    fil=a(indicies);
                    J(2,j)=fil;
                    [X(3,j),X(4,j)]=newMotorPos(Z(1,fil),Z(2,fil),Z(4,fil),Z(5,fil),X(3,j),X(4,j));
                end
                if p0*h > rand() % R leg unattaches
                    J(2,j)=0;
                    X(3,j)=X(1,j);
                    X(4,j)=X(2,j);
                end
            else % L leg can attach if there's a filament in the area
                [count]=findfil2Dv2(Z,X(1,j),X(2,j));
                if p1*h > rand() && sum(count)~=0 % L leg can attach if there's a filament in the area
                    a=find(count);
                    indicies=ceil(rand()*max(size(a)));
                    fil=a(indicies);
                    J(1,j)=fil;
                    [X(1,j),X(2,j)]=newMotorPos(Z(1,fil),Z(2,fil),Z(4,fil),Z(5,fil),X(1,j),X(2,j));
                end
                if p0*h > rand() % L leg unattaches
                    J(1,j)=0;
                    X(1,j)=X(3,j);
                    X(2,j)=X(4,j);
                end
            end
        elseif J(1,j)~=0 && J(2,j) ~=0 && J(1,j)~=J(2,j)
            if rand()>0.5
                if p0*h>rand() % motor will detach
                    J(1,j)=0;
                    X(1,j)=X(3,j);
                    X(2,j)=X(4,j);
                end
            else
                if p0*h>rand()
                    J(2,j)=0;
                    X(3,j)=X(1,j);
                    X(4,j)=X(2,j);
                end
            end
        end
    end
    
    
    % MOTOR MOVEMENT
    for j=1:M
        if J(1,j)==0 && J(2,j)~=0 % L motor is unattached R is attached
            for i=1:N
                if J(2,j)==i % if motor j is attached to filament i
                    S(i,j)=sqrt((Z(1,i)-X(3,j))^2+(Z(2,i)-X(4,j))^2);
                    cb=sqrt((Z(4,i)-X(3,j))^2+(Z(5,i)-X(4,j))^2);
                    if cb>L % then motor has past the end of the filament and will fall off.
                        J(2,j)=0;
                        % place motor at the plus end of the filament
                        X(3,j)=Z(1,i);
                        X(4,j)=Z(2,i);
                        X(1,j)=X(3,j);
                        X(2,j)=X(4,j);
                    else
                        % motor continues to walk down the filament
                        S(i,j)=S(i,j)-v*h;
                        X(3,j)=Z(1,i)-S(i,j)*cos(Z(3,i));
                        X(4,j)=Z(2,i)-S(i,j)*sin(Z(3,i));
                        X(1,j)=X(3,j);
                        X(2,j)=X(4,j);
                    end

                end
            end
            cb=0;
        elseif J(1,j)~=0 && J(2,j)==0 % R motor is unattached L is attached
            for i=1:N
                if J(1,j)==i % if motor j is attached to filament i
                    S(i,j)=sqrt((Z(1,i)-X(1,j))^2+(Z(2,i)-X(2,j))^2);
                    cb=sqrt((Z(4,i)-X(1,j))^2+(Z(5,i)-X(2,j))^2);
                    if cb>L % then motor has past the end of the filament and will fall off.
                        J(1,j)=0;
                        % place motor at the plus end of the filament
                        X(1,j)=Z(1,i);
                        X(2,j)=Z(2,i);
                        X(3,j)=X(1,j);
                        X(4,j)=X(2,j);
                    else
                        % motor continues to walk down the filament
                        S(i,j)=S(i,j)-v*h;
                        X(1,j)=Z(1,i)-S(i,j)*cos(Z(3,i));
                        X(2,j)=Z(2,i)-S(i,j)*sin(Z(3,i));
                        X(3,j)=X(1,j);
                        X(4,j)=X(2,j);
                    end             
                end
            end
            cb=0;
        elseif J(1,j)~=0 && J(2,j)~=0 && J(1,j)~=J(2,j)% R and L motor are attached to two different filaments
            for i=1:N
                if J(1,j)==i % motor j is attached to filament i
                    S(i,j)=sqrt((Z(1,i)-X(1,j))^2+(Z(2,i)-X(2,j))^2);
                    cb=sqrt((Z(4,i)-X(1,j))^2+(Z(5,i)-X(2,j))^2);
                    if cb>L % then motor has past the end of the filament and will fall off.
                        J(1,j)=0;
                        % place motor at the plus end of the filament
                        X(1,j)=Z(1,i);
                        X(2,j)=Z(2,i);
                    else
                        % motor continues to walk down the filament
                        S(i,j)=S(i,j)-v*h;
                        X(1,j)=Z(1,i)-S(i,j)*cos(Z(3,i));
                        X(2,j)=Z(2,i)-S(i,j)*sin(Z(3,i));
                    end
                end
                if J(2,j)==i % motor j is attached to filament i
                    S(i,j)=sqrt((Z(1,i)-X(3,j))^2+(Z(2,i)-X(4,j))^2);
                    cb=sqrt((Z(4,i)-X(3,j))^2+(Z(5,i)-X(4,j))^2);
                    if cb>L % then motor has past the end of the filament and will fall off.
                        J(2,j)=0;
                        % place motor at the plus end of the filament
                        X(3,j)=Z(1,i);
                        X(4,j)=Z(2,i);
                    else
                        % motor continues to walk down the filament
                        S(i,j)=S(i,j)-v*h;
                        X(3,j)=Z(1,i)-S(i,j)*cos(Z(3,i));
                        X(4,j)=Z(2,i)-S(i,j)*sin(Z(3,i));
                    end
                end
            end                    
            motorStretch=sqrt((X(1,j)-X(3,j))^2+(X(2,j)-X(4,j))^2);
            if motorStretch>r % then motor as passed the threshold of stretch and falls off.
                J(1,j)=0;
                J(2,j)=0;
                if rand()>0.5 %equal probability of following off and staying at right as left
                    X(1,j)=X(3,j);
                    X(2,j)=X(4,j);
                else
                    X(3,j)=X(1,j);
                    X(4,j)=X(2,j);
                end
            end
            cb=0;
        elseif J(1,j)==J(2,j) && J(1,j)~=0 %R and L motor are attached to the same filament
            for i=1:N
                if J(1,j)==i %found the filament both sides of motor are attached to
                    S(i,j)=sqrt((Z(1,i)-X(1,j))^2+(Z(2,i)-X(2,j))^2);
                    cb=sqrt((Z(4,i)-X(1,j))^2+(Z(5,i)-X(2,j))^2);
                    if cb>L % then motor has fallen off the filament
                        J(1,j)=0;
                        J(2,j)=0;
                        X(1,j)=Z(1,i);
                        X(2,j)=Z(2,i);
                        X(3,j)=X(1,j);
                        X(4,j)=X(2,j);
                    else
                        S(i,j)=S(i,j)-v*h;
                        X(1,j)=Z(1,i)-S(i,j)*cos(Z(3,i));
                        X(2,j)=Z(2,i)-S(i,j)*sin(Z(3,i));
                        X(3,j)=X(1,j);
                        X(4,j)=X(2,j);
                    end
                end
            end
            cb=0;
        end
%         % CHECK TO SEE IF ANY MOTOR HAS PAST THE BOUNDARY OF THE HEXAGON.
%         % IF SO, DETACH IT AND PUT IT BACK INSIDE THE HEXAGON.
%         inMotor = inpolygon(X(1,j),X(2,j),xhex,yhex);
%         while inMotor==0
%             X(1,j) = xhex(1) + (xhex(4)-xhex(1))*rand();
%             X(2,j) = yhex(2) + (yhex(5)-yhex(2))*rand();
%             inMotor = inpolygon(X(1,j),X(2,j),xhex,yhex);
%             J(1,j)=0;
%             J(2,j)=0;
%         end
%         X(3,j)=X(1,j);
%         X(4,j)=X(2,j);
    end
    
    
    % Filament Polymerization
    for i=1:N %Only occurs for non-vertex filaments
        if p2*h > rand() % If this is true, then depolymerization event occurs causing attached motors to fall off, and filament to change plus end position and angle
            for j=1:M %Find any attached motors to the depolymerizing filament
                if J(1,j)==i && J(2,j) ~=i %if right leg of motor is attached to polymerizing filament i, but left leg is on a different filament
                    J(1,j)=0; % right leg unattaches
                    X(1,j)=X(3,j); %motor moves to the attached left leg position
                    X(2,j)=X(4,j);
                elseif J(2,j)==i  && J(1,j)~=i % if left leg of motor is attached to polymerizing filament i, and right leg is on a different filament
                    J(2,j)=0; % left leg unattaches
                    X(3,j)=X(1,j); % motor moves to the attached right leg position
                    X(4,j)=X(2,j);
                elseif J(1,j)==i && J(2,j)==i %if both legs are on the depolymerizing filament
                    J(1,j)=0; %both legs unattach
                    J(2,j)=0;
                    if rand()>0.5 %equal probability of going to right leg position or left leg position (assuming they could be different in this case)
                        X(1,j)=X(3,j);
                        X(2,j)=X(4,j);
                    else
                        X(3,j)=X(1,j);
                        X(4,j)=X(2,j);
                    end
                end
            end
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
                % has already proven that the plus ends are inside the
                % hexagon, define a new angle of orientation for the
                % filaments, and then redefine a new minus end.
                Z(3,i)=rand()*2*pi;
                Z(4,i)=Z(1,i)-L*cos(Z(3,i));
                Z(5,i)=Z(2,i)-L*sin(Z(3,i));
                inMinus=inpolygon(Z(4,i),Z(5,i),xpol,ypol);
            end
        end
    end
    
    % BOUNDARY POTENTIAL
    % Define a "motor" spring that acts on filaments within a distance
    % BoundaryRadius from the boundary
    vCount=ones(1,max(size(xpol)));
    for i=1:N
        d2=100;
        COMx(i)=(Z(1,i)+Z(4,i))/2;
        COMy(i)=(Z(2,i)+Z(5,i))/2;
        for step=1:max(size(xpol))-1
            d1=(COMx(i)-xpol(step))^2+(COMy(i)-ypol(step))^2;
            if d1<d2
                edgeX(i)=xpol(step);
                edgeY(i)=ypol(step);
                Emod(i)=E(step);
                vCount(step)=max(size(find(edgeX==xpol(step))));
                d2=d1;
            end
        end
    end
    for i=1:N
        rot=[cos(Z(3,i)),sin(Z(3,i)); -sin(Z(3,i)),cos(Z(3,i))];
        for step=1:max(size(xpol))-1
            if xpol(step)==edgeX(i)
                numV=vCount(step);
            end
        end
        if sqrt((edgeX(i)-COMx(i))^2+(edgeY(i)-COMy(i))^2)<BoundaryRadius
            FBx(i)=Emod(i)*th*perim*(edgeX(i)-Z(1,i))/(Np*numV);
            FBy(i)=Emod(i)*th*perim*(edgeY(i)-Z(2,i))/(Np*numV);
        else
            FBx(i)=0;
            FBy(i)=0;
            angFB(i)=0;
        end
        FB(1,i)=rot(1,1)*FBx(i)+rot(1,2)*FBy(i);
        FB(2,i)=rot(2,1)*FBx(i)+rot(2,2)*FBy(i);
        angFB(i)=sqrt((edgeX(i)-Z(1,i))^2+(edgeY(i)-Z(2,i))^2)*FB(2,i);
    end
            
    
    
    % FILAMENT MOVEMENT   
    % movement for free filaments due to bound motors
    for i=1:N
        rot=[cos(Z(3,i)),sin(Z(3,i)); -sin(Z(3,i)),cos(Z(3,i))];
        rot_inv=[cos(Z(3,i)),-sin(Z(3,i));sin(Z(3,i)),cos(Z(3,i))];
%         COMx=(Z(1,i)+Z(4,i))/2;
%         COMy=(Z(2,i)+Z(5,i))/2;
        XR=rot*[COMx(i);COMy(i)];
        len=zeros(1,M);
        fx=zeros(1,M);
        fy=zeros(1,M);
        for j=1:M
            if J(1,j)==i && J(1,j)~=J(2,j) && J(2,j) ~=0 % both ends of motor are attached to 2 different filaments
                mxA(j)=X(1,j);
                myA(j)=X(2,j);
                mxB(j)=X(3,j);
                myB(j)=X(4,j);
                len(j)=sqrt((mxA(j)-COMx(i))^2+(myA(j)-COMy(i))^2);
                if S(i,j)>L/2;
                    len(j)=-len(j);
                end
                x1(j)=COMx(i)+len(j)*cos(Z(3,i));
                y1(j)=COMy(i)+len(j)*sin(Z(3,i));
                fx(j)=k*(mxB(j)-x1(j));
                fy(j)=k*(myB(j)-y1(j));
            elseif J(2,j)==i && J(1,j)~=J(2,j) && J(1,j)~=0
                mxA(j)=X(3,j);
                myA(j)=X(4,j);
                mxB(j)=X(1,j);
                myB(j)=X(2,j);
                len(j)=sqrt((mxA(j)-COMx(i))^2+(myA(j)-COMy(i))^2);
                if S(i,j)>L/2;
                    len(j)=-len(j);
                end
                x1(j)=COMx(i)+len(j)*cos(Z(3,i));
                y1(j)=COMy(i)+len(j)*sin(Z(3,i));
                fx(j)=k*(mxB(j)-x1(j));
                fy(j)=k*(myB(j)-y1(j));
            end
            F(1,j)=fx(j);
            F(2,j)=fy(j);
%             FR(:,j)=rot*F(:,j);
            FR(1,j)=rot(1,1)*F(1,j)+rot(1,2)*F(2,j);
            FR(2,j)=rot(2,1)*F(1,j)+rot(2,2)*F(2,j);
            ang(j)=len(j)*FR(2,j);
        end

        
        %Update positions
        XR_n(1,1)=XR(1)+h/l_pr*(sum(FR(1,:))+FB(1,i));
        XR_n(2,1)=XR(2)+h/l_pp*(sum(FR(2,:))+FB(2,i));
        alpha_n=Z(3,i)+h/l_r*(sum(ang)+angFB(i));

        %transform back
        COM_n=rot_inv*XR_n;
        Z(3,i)=alpha_n;
        Z(1,i)=COM_n(1,1)+(L/2)*cos(Z(3,i));
        Z(2,i)=COM_n(2,1)+(L/2)*sin(Z(3,i));
        Z(4,i)=Z(1,i)-L*cos(Z(3,i));
        Z(5,i)=Z(2,i)-L*sin(Z(3,i));

        %update motors to stay on the newly relocated filament
        for j=1:M
            if J(1,j)~=J(2,j) && J(1,j)~=0 && J(2,j)~=0 % both ends of motor are attached to different filaments
                if J(1,j)==i 
                    X(1,j)=Z(1,i)-S(i,j)*cos(Z(3,i));
                    X(2,j)=Z(2,i)-S(i,j)*sin(Z(3,i));
                elseif J(2,j)==i
                    X(3,j)=Z(1,i)-S(i,j)*cos(Z(3,i));
                    X(4,j)=Z(2,i)-S(i,j)*sin(Z(3,i));
                end
            end
        end
%         clear rot rot_inv COMx COMy XR mxA myA mxB myB F FR P_n COM_n
    end
                    
    
     

    cd(newFolder);
    fid=fopen(sprintf('fil%d.txt',t),'w');
    for i=1:N
        fprintf(fid,'%f  %f  %f\n',[Z(1,i),Z(2,i),Z(3,i)]);
    end
    fclose(fid);
    clear fid;
    
    fidMot=fopen(sprintf('mot%d.txt',t),'w');
    for j=1:M
        fprintf(fidMot,'%f  %f  %f  %f  %f  %f\n',[X(1,j),X(2,j),X(3,j),X(4,j),J(1,j),J(2,j)]);
    end
    fclose(fidMot);
    clear fidMot;
    cd(runSims);


end


