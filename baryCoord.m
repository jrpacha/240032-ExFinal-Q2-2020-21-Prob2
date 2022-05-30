function [alphas, isInside]=baryCoord(vertices, point)
% (c)Numerical Factory 2021
isInside=1; 
A=[ones(3,1), vertices];
b=[1;0;0];
c1=A\b; 
b=[0;1;0];
c2=A\b;
b=[0;0;1];
c3=A\b;

Psi1 = @(x,y) c1(1)+c1(2)*x+c1(3)*y;
Psi2 = @(x,y) c2(1)+c2(2)*x+c2(3)*y;
Psi3 = @(x,y) c3(1)+c3(2)*x+c3(3)*y;
%
% Compute barycentric coordinates
%
alpha1 = Psi1(point(1),point(2));
alpha2 = Psi2(point(1),point(2));
alpha3 = Psi3(point(1),point(2));
alphas=[alpha1,alpha2,alpha3];
if ( min(alphas) < -1.e-14) % add max(alphas) > 1 is not needed
    isInside=0;
end