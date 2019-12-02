%RoseDiagramForce.m
% 11/18/14
% Callie J Miller
%
% The purpose of this code is to create a rose diagram of the boundary
% filament's forces at specific time points.
clear all
close all
old=cd;
%cd('/Users/CallieMiller/Desktop/Chanet, Martin, and Us Manuscript Data/Simulations/100-gon/Soft Top Stiff Sides');
% cd('/Users/CallieMiller/Desktop/Chanet, Martin, and Us Manuscript Data/Simulations/100-gon/Stiff Top Soft Sides');
%cd('/Users/CallieMiller/Desktop/Chanet, Martin, and Us Manuscript Data/Simulations/100-gon/Control, 1 kPa');
% cd('/Users/CallieMiller/Desktop/Chanet, Martin, and Us Manuscript Data/Simulations/100-gon/SoftTop (50 Pa) Stiff Sides');
folder_name=uigetdir;
cd(folder_name);


% Defining vertex parameters
E=zeros(1,100);
for n=1:100
    if (n<15 || n>89 || (n>39 && n<65)) % Sides of 100 sided polygon
        E(n)=20; 
    else
        E(n)=100; 
    end
end

N=1000;
hrad=2; %radius of patch
vrad=2;
x0=0; %x component of center
y0=0; %y component of center
cd(old);
[xpol,ypol,POLY]=ellipse(hrad,vrad,x0,y0,100);

L=1;
BoundaryRadius=0.5;
th=2;
perim=2*pi*hrad;


for t=1:1000
Force=zeros(3,100);
Xvertex=zeros(1,100);
Yvertex=zeros(1,100);
count=ones(1,100);
cd(folder_name);


    empty=zeros(1,100);
    %Read in plus-end positions
    fid=fopen(sprintf('fil%d.txt',t));
    A=fscanf(fid,'%g',[3,inf]);
    fclose(fid);
    
    for i=1:N
        Z(1,i)=A(1,i);
        Z(2,i)=A(2,i);
        Z(3,i)=A(3,i);
        Z(4,i)=Z(1,i)-L*cos(Z(3,i));
        Z(5,i)=Z(2,i)-L*sin(Z(3,i));
    end
    
    clear A fid
    
    % Determine which vertex point is closeest to COM of filaments
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
    
    %If distance between COM and edge is less than BoundaryRadius, then add
    %it to sum up all "springs" at a vertex point acting on that vertex
    for i=1:N
        if sqrt((edgeX(i)-COMx(i))^2+(edgeY(i)-COMy(i))^2)<BoundaryRadius
            a=find(edgeX(i)==xpol,max(size(xpol))-1);
            b=find(edgeY(i)==ypol,max(size(xpol))-1);
            ind=a(find(a==b));
            Xvertex(ind)=Xvertex(ind)+Z(1,i);
            Yvertex(ind)=Yvertex(ind)+Z(2,i);
            count(ind)=count(ind)+1;
            empty(ind)=1;
        end
    end
    
    % Calculate the force on each vertex point
    for v=1:100
        Force(2,v)=empty(v)*E(v)*th*perim*(sqrt((Xvertex(v)-xpol(v))^2+(Yvertex(v)-ypol(v))^2))/(100*(count(v)));
    end
      
    Force(3,1:10)=mean(Force(2,1:9));
    Force(3,11:20)=mean(Force(2,11:20));
    Force(3,21:30)=mean(Force(2,21:30));
    Force(3,31:40)=mean(Force(2,31:40));
    Force(3,41:50)=mean(Force(2,41:50));
    Force(3,51:60)=mean(Force(2,51:60));
    Force(3,61:70)=mean(Force(2,61:70));
    Force(3,71:80)=mean(Force(2,71:80));
    Force(3,81:90)=mean(Force(2,81:90));
    Force(3,91:100)=mean(Force(2,91:100));
    
    fid=fopen(sprintf('Rose%d.txt',t),'w');
    for i=1:100
        fprintf(fid,'%f %f %f\n',[Force(1,i),Force(2,i),Force(3,i)]);
    end
    fclose(fid);
    clear fid Z Force COMx COMy empty edgeX edgeY Xvertex Yvertex
end