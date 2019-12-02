function []=VariableTestingSimMulti(MainrunSims,MainnewFolder,Stiffness)
% VariableTestingSim.m
% 01/26/14
% Callie J Miller
% 
% The purpose of this code is to allow the user to quickly test various
% variables of the simulation of actomyosin in a hexagon.
% USES THE FOLLOWING FUNCTIONS: initialize.m, HexSim.m (in the Simulation
% Codes folder)
% MODIFIED 03/02/14: Test if can get angular co-alignment by manipulating
% vertex filament spring stiffnesses or lengths.
%
% UPDATED 05/01/15: To account for incorrectness with force calculations
%
% UPDATED 02/12/16: Versus boundary filaments, implement a potential force
% at the boundary.

global runSims newFolder th BoundaryRadius Time M N L p0 p1 p2 h r v k Np cote x0 y0 xhex yhex Zold distSide maxX pointxm maxY pointyM minY pointym

runSims=MainrunSims;
newFolder=MainnewFolder;

% TimePrompt = 'How long do you want the simulation to run?';
% Time = input(TimePrompt);
% 
% MPrompt='How many motors?';
% M=input(MPrompt);
% 
% NPrompt='How many filaments?';
% N=input(NPrompt);

Time=1000;
M=5000;
N=1000;




%Variables already predeclared
%Time=200;
%M=100;
%N=50;

L=1; %length of filaments
p0=1; %Detachment rate
p1=10; %Attachment rate
p2=0.7;%0.7; %Rate of polymerization 
h=0.01; %stepsize
r=0.3; %radius motor's will attach
% r=2/10; %radius of searching for motor to attach and maximum stretch of motor allowed (the hexagonal radius divided by 10)
v=1; %speed of motor
k=3; %spring stiffness of motor

th=2; % thickness of cortex approximately 2 um
BoundaryRadius=0.5; %assume boundary potential acts on filaments within 0.5 um from boundary edge


Np=100; % number of vertex points to define boundary

tic
[J,Z,X]=initialize();
cd(newFolder);
fid=fopen('fil0.txt','w');
    for i=1:N 
        fprintf(fid,'%f %f %f\n',[Z(1,i),Z(2,i),Z(3,i)]);
    end
    fclose(fid);
    clear fid
    
    fidMot=fopen('mot0.txt','w');
    for j=1:M
            fprintf(fidMot,'%f  %f  %f  %f  %f  %f\n',[X(1,j),X(2,j),X(3,j),X(4,j),J(1,j),J(2,j)]);
    end
    fclose(fidMot);
    cd(runSims);
[X,J,Z,S]=HexSimPolyMulti(X,J,Z,Stiffness);
toc