function prob3e(ver_ponto)

warning('off','all')
%% min(x1^4-2*x1^2*x2+x1^2+x1*x2^2-2*x1+4)
% s.a.: .25*x1^2+.75*x2^2-1<=0; 2*x1^2+x2^2=0

syms lambdaz 
syms x1 x2 real

thetas=x1^4-2*x1^2*x2+x1^2+x1*x2^2-2*x1+4; %funcao a minimizar
thetas_desig=.25*x1^2+.75*x2^2-1;ndes=size(thetas_desig,1); %sujeito a desigualdade
thetas_ig=2*x1^2+x2^2-2;ni=size(thetas_ig,1); %sujeito a igualdade

gradfx=sym(zeros(2,1));

% Gradiente
dfdx1=diff(thetas,x1);dfdx2=diff(thetas,x2);
gradfx(1,1)=dfdx1;gradfx(2,1)=dfdx2;

%gradientes das igualdades/desigualdades

gradg=sym(zeros(2,ndes));
for j=1:ndes
gradg(:,j)=[diff(thetas_desig(j),x1);diff(thetas_desig(j),x2)];
end

gradi=sym(zeros(2,ni));
for j=1:ni
gradi(:,j)=[diff(thetas_ig(j),x1);diff(thetas_ig(j),x2)];
end

%escolher os multiplicadores de lagrange do problema
syms u1 v1 u2 v2 u3 v3

u=[u1]; nu=size(u,2); v=[v1]; nv=size(v,1);

sumg=0;sumg1=0;
    for i=1:nu
    sumg=sumg+u(i)*gradg(:,i)';
    sumg1=sumg1+u(i)*thetas_desig(:,i)';
    end
    somaG=vpa(sumg);
    sumg1=vpa(sumg1);
    
sumh=0;
    for i=1:nv
    sumg=sumh+v(i)*gradi(:,i)';
    end
    somaH=vpa(sumg);

eq(1,1)=vpa(gradfx(1,1)+somaG(1)+somaH(1));
eq(2,1)=vpa(gradfx(2,1)+somaG(2)+somaH(2));

%substituir nos gradientes para verificar as condicoes em pontos especificos

gf(1,1)=double(subs(gradfx(1,1),[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]));gf(1,2)=double(subs(gradfx(2,1),[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]));

for ii=1:ndes
gg(1,ii)=double(subs(gradg(1,ii),[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]));gg(2,ii)=double(subs(gradg(2,ii),[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]));
end

for ii=1:ni
gh(1,ii)=double(subs(gradi(1,ii),[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]));gh(2,ii)=double(subs(gradi(1,ii),[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]));
end

eq1=subs(eq(1,1),[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]);
eq2=subs(eq(2,1),[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]);

[su1,sv1]=solve(eq1,eq2,u1,v1);
if ~strcmp(class(su1),'sym');su1=double(su1);end
if ~strcmp(class(su1),'sym');sv1=double(sv1);end

kktcnd3=subs(sumg1,[x1,x2,u1],[ver_ponto(1,1),ver_ponto(1,2),su1]);
if ~strcmp(class(kktcnd3),'sym');kktcnd3=double(kktcnd3);end

if (su1>=0 & ~isempty(su1)) & ~isempty(sv1) & kktcnd3==0
    fprintf('-----------------------------------------------------------------------------\n')
    fprintf('O ponto (%6.3f,%6.3f) satisfaz as condicoes de KKT!\n',ver_ponto(1,1),ver_ponto(1,2));
    fprintf('u1=%6.3f;v1=%6.3f\n',su1,sv1);
    fprintf('-----------------------------------------------------------------------------\n')
else
    fprintf('-----------------------------------------------------------------------------\n')
    fprintf('O ponto (%6.3f,%6.3f) nao satisfaz as condicoes de KKT!\n',ver_ponto(1,1),ver_ponto(1,2));
    fprintf('-----------------------------------------------------------------------------\n')
end

end