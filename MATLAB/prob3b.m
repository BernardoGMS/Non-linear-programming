function [x,miu,miv,v,u,thetaz,igthetaz,desigthetaz,alpha]=prob3b

%% Formulacao Langrangeano Aumentado

syms x1 x2 real
mf=1;
if mf==1 %Exercicio proposto 3
thetaz=x1^4-2*x1^2*x2+x1^2+x1*x2^2-2*x1+4; %funcao a minimizar
desigthetaz=.25*x1^2+.75*x2^2-1;ndes=size(desigthetaz,1); %sujeito a desigualdade
igthetaz=2*x1^2+x2^2-2;nig=size(igthetaz,1); %sujeito a igualdade
elseif mf==2 %Exemplo do livro
thetaz=(x1-6)^2+(x2-7)^2; %funcao a minimizar
desigthetaz=[-3*x1-2*x2+6;-x1+x2-3;x1+x2-7;2/3*x1-x2-4/3];ndes=size(desigthetaz,1); %sujeito a desigualdade
igthetaz=0;nig=size(igthetaz,1); %sujeito a igualdade
elseif mf==3 %Exercicio proposto 4
thetaz=x1^2-x1*x2+2*x2^2-4*x1-5*x2; %funcao a minimizar
desigthetaz=[x2-2;x1+2*x2-6];ndes=size(desigthetaz,1); %sujeito a desigualdade
igthetaz=0;nig=size(igthetaz,1); %sujeito a igualdade
elseif mf==4 %Exemplo do livro
thetaz=x1^2+x2^2; %funcao a minimizar
desigthetaz=[-x1-x2+4;-2*x1-x2+5];ndes=size(desigthetaz,1); %sujeito a desigualdade
igthetaz=0;nig=size(igthetaz,1); %sujeito a igualdade
end

syms rhou rhov v1 u1
%% Passo de initicializacao
miu=ones(1,ndes)*.01;miv=ones(1,nig)*10;
v=zeros(1,nig);u=zeros(1,ndes);

alpha=thetaz+miu*(desigthetaz+u'./(2*miu')).^2-sum(u.^2./(4*miu))+v*igthetaz+miv*igthetaz.^2;
alphass=thetaz+rhou*(desigthetaz+u1/(2*rhou)).^2-u1^2./(4*rhou)+v1*igthetaz+rhov*igthetaz.^2;
gradiente=[diff(alpha,x1) diff(alpha,x2)];
fprintf('Formulacao de Lagrangeano Aumentado:\n')
disp(alphass)
[xs1,xs2]=solve(gradiente(1),gradiente(2));
x(1,:)=[double(vpa(xs1(1,1))) double(vpa(xs2(1,1)))];

%% Preparar o plot da funcao
xp1=-1.5:0.01:2.5;xp2=-1.5:0.01:2.5;
[X1,X2]=meshgrid(xp1,xp2);

if mf==1
z=X1.^4-2*X1.^2*X2+X1.^2+X1.*X2.^2-2*X1+4;
dz1=.25*X1.^2+.75*X2.^2-1;
iz1=2*X1.^2+X2.^2-2;
elseif mf==2
z=(X1-6).^2+(X2-7).^2;
dz1=-3*X1-2*X2+6;
dz2=-X1+X2-3;
dz3=X1+X2-7;
dz4=2/3*X1-X2-4/3;
elseif mf==3
z=X1.^2-X1.*X2+2*X2.^2-4*X1-5*X2;
dz1=X2-2;dz2=X1+2*X2-6; 
elseif mf==4
z=X1.^2+X2.^2;
dz1=-X1-X2+4;
dz2=-2*X1-X2+5;
end

titulo='Minimizacao - Lagrange';
            
scrsz = get(0,'ScreenSize');
posfig=[scrsz(3) scrsz(4) scrsz(3) scrsz(4)];

figure('Name',['Problema 2: ',titulo],...
'NumberTitle','off','OuterPosition',posfig);
hold on;

line([0 2.5],[0 0],'Color','k', ...
    'LineWidth',2,'LineStyle','-');
line([0 0],[0 2.5],'Color','k', ...
    'LineWidth',2,'LineStyle','-');

if mf==1
contour(X1,X2,z,'ShowText','on');hold on;
contour(X1,X2,dz1,[0 0],'m-');contour(X1,X2,iz1,[0 0],'g-');
xlabel('X1');ylabel('X2');hold on;legend('x1>=0','x2>=0','f=x_1^4-2*x_1^2*x_2+x_1^2+x_1*x_2^2-2*x_1+4','0>=.25*x_1^2+.75*x_2^2-1','2*x_1^2+x_2^2-2=0')
elseif mf==2
contour(X1,X2,z,'ShowText','on');hold on;
contour(X1,X2,dz1,[0 0],'r-');contour(X1,X2,dz2,[0 0],'r-');contour(X1,X2,dz3,[0 0],'r-');contour(X1,X2,dz4,[0 0],'r-');
xlabel('X1');ylabel('X2');hold on;
elseif mf==3
contour(X1,X2,z,'ShowText','on');hold on;
contour(X1,X2,dz1,[0 0],'r-');contour(X1,X2,dz2,[0 0],'r-');
xlabel('X1');ylabel('X2');hold on;
elseif mf==4
contour(X1,X2,z,'ShowText','on');hold on;
contour(X1,X2,dz1,[0 0],'r-');contour(X1,X2,dz2,[0 0],'r-');
xlabel('X1');ylabel('X2');hold on;
end
grid on
    
end