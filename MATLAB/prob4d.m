function [x,miu,miv,thetaz,igthetaz,desigthetaz,alpha]=prob4d

    mf=3;%Exercicio proposto 4
syms x1 x2 real
if mf==1 %Exercicio proposto 3
thetaz=x1^4-2*x1^2*x2+x1^2+x1*x2^2-2*x1+4; %funcao a minimizar
desigthetaz=.25*x1^2+.75*x2^2-1;ndes=size(desigthetaz,1); %sujeito a desigualdade
igthetaz=2*x1^2+x2^2-2;ni=size(igthetaz,1); %sujeito a igualdade
elseif mf==2 %Exemplo do livro
thetaz=(x1-6)^2+(x2-7)^2; %funcao a minimizar
desigthetaz=[-3*x1-2*x2+6;-x1+x2-3;x1+x2-7;2/3*x1-x2-4/3];ndes=size(desigthetaz,1); %sujeito a desigualdade
igthetaz=0;ni=size(igthetaz,1); %sujeito a igualdade
elseif mf==3 %Exercicio proposto 4
thetaz=x1^2-x1*x2+2*x2^2-4*x1-5*x2; %funcao a minimizar
desigthetaz=[x1+2*x2-6;x2-2];ndes=size(desigthetaz,1); %sujeito a desigualdade
igthetaz=0;ni=size(igthetaz,1); %sujeito a igualdade
elseif mf==4 %Exemplo do livro
thetaz=x1^2+x2^2; %funcao a minimizar
desigthetaz=[-x1-x2+4;-2*x1-x2+5];ndes=size(desigthetaz,1); %sujeito a desigualdade
igthetaz=0;ni=size(igthetaz,1); %sujeito a igualdade
end

syms rhou1 rhou2
%% Passo de initicializacao
miu=zeros(1,ndes)/10;miv=zeros(1,ni)/10;

alpha=thetaz+miu*desigthetaz.^2+miv*igthetaz.^2;
alphass=thetaz+[rhou1 rhou2]*(desigthetaz.^2); %apenas para o problema 4
fprintf('Formulacao de Penalidade:\n')
disp(alphass);fprintf('Sendo rhou1 e rhou2 as penalidades dos contrangimentos 1 e 2, respetivamente\n')
gradiente=[diff(alpha,x1) diff(alpha,x2)];
[xs1,xs2]=solve(gradiente(1),gradiente(2));
x(1,:)=[double(vpa(xs1(1,1))) double(vpa(xs2(1,1)))];

%% Preparar o plot da funcao

if mf==1
xp1=0:0.01:3;xp2=0:0.01:3;
[X1,X2]=meshgrid(xp1,xp2);
z=X1.^4-2*X1.^2*X2+X1.^2+X1.*X2.^2-2*X1+4;
dz1=.25*X1.^2+.75*X2.^2-1;
iz1=2*X1.^2+X2.^2-2;
elseif mf==2
xp1=0:0.01:8;xp2=0:0.01:8;
[X1,X2]=meshgrid(xp1,xp2);
z=(X1-6).^2+(X2-7).^2;
[x1dz1,y1dz1]=solve(desigthetaz(1,1),desigthetaz(2,1),x1,x2);[x3dz1,y3dz1]=solve(desigthetaz(1,1),desigthetaz(4,1),x1,x2);
[x1dz2,y1dz2]=solve(desigthetaz(2,1),desigthetaz(3,1),x1,x2);[x1dz3,y1dz3]=solve(desigthetaz(3,1),desigthetaz(4,1),x1,x2);
xpp=double([x1dz1,x3dz1,x1dz2,x1dz3]);ypp=double([y1dz1,y3dz1,y1dz2,y1dz3]);
x1max=max(xpp);x2max=max(ypp);xp1=0:0.01:x1max;xp2=0:0.01:x2max;[X1i,X2i]=meshgrid(xp1,xp2);
dz1=-3*X1i-2*X2i+6;
dz2=-X1i+X2i-3;
dz3=X1i+X2i-7;
dz4=2/3*X1i-X2i-4/3;
elseif mf==3
xp1=0:0.01:3.5;xp2=0:0.01:3.5;
[X1,X2]=meshgrid(xp1,xp2);
z=X1.^2-X1.*X2+2*X2.^2-4*X1-5*X2;
[x1dz1,y1dz1]=solve(desigthetaz(1,1),desigthetaz(2,1),x1,x2);
xpp=double(x1dz1);ypp=double(y1dz1);
x1max=max(xpp);x2max=max(ypp);xp1=0:0.01:x1max;xp2=0:0.01:x2max;[X1i,X2i]=meshgrid(xp1,xp2);
dz1=X2i-2;x1max=max(xpp);xp1=x1max:0.01:3.5;xp2=0:0.01:3.5;[X1ii,X2ii]=meshgrid(xp1,xp2);
dz2=X1ii+2*X2ii-6; 
elseif mf==4
xp1=0:0.01:5;xp2=0:0.01:5;
[X1,X2]=meshgrid(xp1,xp2);
z=X1.^2+X2.^2;
dz1=-X1-X2+4;
dz2=-2*X1-X2+5;
end

titulo='Minimizacao - Formulacao de Penalidade';
            
scrsz = get(0,'ScreenSize');
posfig=[scrsz(3) scrsz(4) scrsz(3) scrsz(4)];

figure('Name',['Problema 4: ',titulo],...
'NumberTitle','off','OuterPosition',posfig);
hold on;

if mf==1
contour(X1,X2,z,'ShowText','on');hold on;
contour(X1,X2,dz1,[0 0],'r-','Linewidth',4);contour(X1,X2,iz1,[0 0],'g-','Linewidth',4);
xlabel('X1');ylabel('X2');hold on;
elseif mf==2
contour(X1,X2,z,'ShowText','on');hold on;
contour(X1i,X2i,dz1,[0 0],'r-','Linewidth',4);contour(X1i,X2i,dz2,[0 0],'r-','Linewidth',4);contour(X1i,X2i,dz3,[0 0],'r-','Linewidth',4);contour(X1i,X2i,dz4,[0 0],'r-','Linewidth',4);
xlabel('X1');ylabel('X2');hold on;
elseif mf==3
contour(X1,X2,z,'ShowText','on');hold on;[C1,h1] = contour(X1,X2,z,-10.78,'g--');clabel(C1,h1);
contour(X1i,X2i,dz1,[0 0],'b-','Linewidth',4);contour(X1ii,X2ii,dz2,[0 0],'r-','Linewidth',4);
xlabel('X1');ylabel('X2');hold on;
line([0 3],[0 0],'Color','k','LineWidth',4,'LineStyle','-');line([0 0],[0 2],'Color','k','LineWidth',4,'LineStyle','-');
legend('f=x_1^2-*x_1*x_2+2*x_2-4x_1-5*x_2','f(x_o_p_t)','0>=x_2-2','0>=x_1+2x_2-6','x1>=0','x2>=0')
elseif mf==4
contour(X1,X2,z,'ShowText','on');hold on;
contour(X1,X2,dz1,[0 0],'r-');contour(X1,X2,dz2,[0 0],'r-');
xlabel('X1');ylabel('X2');hold on;
end

end