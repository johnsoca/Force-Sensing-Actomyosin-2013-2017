function [T9]=findfil2Dv2(Z,x1,y1)
global N r

D=zeros(1,N);
D1=zeros(1,N);
Dpp=zeros(1,N);
Dmm=zeros(1,N);
% Logical Method:
for i=1:N
    D(i)=((x1-Z(1,i))*(Z(4,i)-Z(1,i))+(y1-Z(2,i))*(Z(5,i)-Z(2,i)))/((Z(4,i)-Z(1,i))^2+(Z(5,i)-Z(2,i))^2);
end
T1=logical(D>0);
TT=logical(D<1);
T1=logical(T1+TT==2);
for i=1:N
    if T1(i)==1;
        Cx=x1-Z(1,i);
        Cy=y1-Z(2,i);
        Dx=Z(4,i)-Z(1,i);
        Dy=Z(5,i)-Z(2,i);
        D1(i)=(Cx*Dy-Dx*Cy)^2/(Dx^2+Dy^2);
    else
        Dpp(i)=(x1-Z(1,i))^2+(y1-Z(2,i))^2;
        Dmm(i)=(x1-Z(4,i))^2+(y1-Z(5,i))^2;
    end
end
T3=logical(D1~=0);
T2=logical(D1<r^2);
T4=logical(T3+T2==2);

T5=logical(Dpp<r^2);
T6=logical(Dmm<r^2);

T7=T5+T6;
T8=logical(T7==1);
T9=T4+T8; %when T4 is true and when T5 or T6 are true)