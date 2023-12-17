function prob4b

warning('off','all')

syms lambdaz 
syms x1 x2 real

thetas=x1^2-x1*x2+2*x2^2-4*x1-5*x2; %funcao a minimizar
thetas_desig=[x1+2*x2-6;x2-2];ndes=size(thetas_desig,1); %sujeito a desigualdade
thetas_ig=0;ni=size(thetas_ig,1); %sujeito a igualdade

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

u=[u1,u2]; nu=size(u,2); v=[v1]; nv=size(v,2);

sumg=0;sumg1=sym(zeros(nu,1));
    for i=1:nu
    sumg=sumg+u(i)*gradg(:,i)';
    sumg1(i,1)=u(i)*thetas_desig(i,1);
    end
    somaG=vpa(sumg);
    sumg1=vpa(sumg1);
    
sumh=0;
    for i=1:nv
    sumg=sumh+v(i)*gradi(:,i)';
    end
    somaH=vpa(sumg);
    
eq1=gradfx(1,1)+somaG(1)+somaH(1);
eq2=gradfx(2,1)+somaG(2)+somaH(2);

fprintf('\n')
fprintf('-----------------------------------------------------------------------------\n')
fprintf('Condicoes necessarias de KKT:\n%s=0\n%s=0\n',char(eq1),char(eq2));
fprintf('%s=%g\n%s=%g\n',char(sumg1(1,1)),0,char(sumg1(2,1)),0);
for ii=1:nu
fprintf('%s>=%g\n',char(u(ii)),0);
end
fprintf('-----------------------------------------------------------------------------\n')
fprintf('\n')


end