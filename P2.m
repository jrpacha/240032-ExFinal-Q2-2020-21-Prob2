clearvars;
close all
clc

%Data
kc = 1.1;
beta = 21;
xp = 3.6;

tempInner = 100.0;
tempInf = 30;


%Geometria
eval('chimenea2');
numNod=size(nodes,1);
numElem=size(elem,1);
numbering=0;


%Find nodes in the outer boundary
nodsL = find(nodes(:,1) < 0.001);
nodsR = find(nodes(:,1) > 19.999);
nodsB = find(nodes(:,2) < 0.001);
nodsT = find(nodes(:,2) > 19.999);

nodsOuter = unique([nodsL; nodsR; nodsB; nodsT]);

%Find the nodes in the inner boundary
nodsInnerL = find(abs(nodes(:,1)-6) < 0.001 & abs(nodes(:,2)-10) < 4.001);
nodsInnerR = find(abs(nodes(:,1)-14) < 0.001 & abs(nodes(:,2)-10) < 4.001);
nodsInnerB = find(abs(nodes(:,1)-10) < 4.001 & abs(nodes(:,2)-6) < 0.001);
nodsInnerT = find(abs(nodes(:,1)-10) < 4.001 & abs(nodes(:,2)-14) < 0.001);

nodsInner = unique([nodsInnerL; nodsInnerR; nodsInnerB; nodsInnerT]); 
plotElementsOld(nodes,elem,numbering);
hold on
plot(nodes(nodsInner,1),nodes(nodsInner,2),'o','markerFaceColor','red','markerSize',10);
plot(nodes(nodsOuter,1),nodes(nodsOuter,2),'o','markerFaceColor','green','markerSize',10);
hold off

%Coeficients
a11=kc;
a12=0;
a21=a12;
a22=a11;
a00=0;
f=0;
coeff=[a11,a12,a21,a22,a00,f];

%Global system
K=zeros(numNod);   %global Stiff Matrix
F=zeros(numNod,1); %global internal forces vector
Q=zeros(numNod,1); %global secondary variables vector

for e=1:numElem 
  [Ke,Fe]=linearTriangElement(coeff,nodes,elem,e);
  %
  % Assemble the elements
  %
  rows=[elem(e,1); elem(e,2); elem(e,3)];
  colums= rows;     
  K(rows,colums)=K(rows,colums)+Ke; %assembly
  if (coeff(6) ~= 0)
      F(rows)=F(rows)+Fe;
  end
end % end for elements
% we save a copy of the initial F array
% for the postprocess step
Kini=K;
Fini=F;

fixedNods = nodsInner;
freeNods = setdiff(1:numNod,fixedNods);

%Boundary Conditions
%Natural (convection)        
[K,Q]=applyConvTriang(nodsOuter',beta,tempInf,K,Q,nodes,elem);
% [K,Q]=applyConvTriang(nodsL',beta,tempInf,K,Q,nodes,elem);
% [K,Q]=applyConvTriang(nodsR',beta,tempInf,K,Q,nodes,elem);
% [K,Q]=applyConvTriang(nodsB',beta,tempInf,K,Q,nodes,elem);
% [K,Q]=applyConvTriang(nodsT',beta,tempInf,K,Q,nodes,elem);


%Essential
u = zeros(numNod,1);
u(nodsInner) = tempInner;

%Reduced system
Qm = F(freeNods) + Q(freeNods) - K(freeNods,fixedNods)*u(fixedNods);
Km = K(freeNods,freeNods);

%Solve the reduced system
um = Km\Qm;

u(freeNods) = um;

%Post process
QAfterPostP = Kini*u-Fini;

%Contour plot
title = 'Contour plot of the temperature';
colorScale = 'jet';
plotContourSolution(nodes,elem,u,title,colorScale)

q = [xp, 10;
    20-xp,10;
    10,xp;
    10,20-xp];

interpTemp = zeros(4,1); %To hold the interpolated temperature at the
                         %four given points

for s = 1:4
    point = q(s,:);
    for e=1:numElem
        nodesElem = elem(e,:);
        vertexs= nodes(nodesElem,:);
        [alphas,isInside] = baryCoord(vertexs,point);
        if (isInside > 0)
            interpTemp(s) = alphas*u(nodesElem);
            break;
        end
    end
end

%Write the solutions
fprintf('(a) K(1,1) before applying the convective B.C: K(1,1) = %.5e\n',...
    Kini(1,1))
fprintf('    Check. K(5,5) = %.5e (before applying the convective BC).\n',...
    Kini(5,5))
fprintf('(b) Q(1) after post-process, Q(1) = %.5e\n',QAfterPostP(1));
fprintf('   Check. Its value before the post-process is Q(1) = %.5e\n',...
    Q(1))
fprintf('(c) max(u(j)) - min(u(j)) = %.5e\n',...
    max(interpTemp)-min(interpTemp));
fprintf('    Check. The interpolated temperature u(1) at q(1) is %.5e\n',...
    interpTemp(1))