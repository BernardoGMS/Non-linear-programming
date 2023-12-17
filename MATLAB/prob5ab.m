function thetas=prob5ab

%% min((2*x^2+x*y+y^2+y*z+x^2-6*x-7*y-8*z+9)

syms x y z ZZ real
syms lambdaz
thetas=2*x^2+x*y+y^2+y*z+z^2-6*x-7*y-8*z+9; %funcao a minimizar

%% a) Condicoes necessarias de minimo:

% Condicao de 1ordem para x* ser um optimo local:
% E d tal que grad(f(x*))t.d<0 --> E delta>0 tal que f(x*+lambda.d)<f(x*)
%para cada lambda pertencente [0, delta];

gradfx=sym(zeros(3,1));hessian=sym(zeros(3,3));
dfdx=diff(thetas,x);dfdy=diff(thetas,y);dfdz=diff(thetas,z);
df2dx=diff(dfdx,x);df2dy=diff(dfdy,y);df2dz=diff(dfdz,z);
df2dxy=diff(dfdy,x);df2dxz=diff(dfdz,x);df2dyz=diff(dfdz,y);
gradfx(1,1)=dfdx;gradfx(2,1)=dfdy;gradfx(3,1)=dfdz;
hessian(1,1)=df2dx;hessian(1,2)=df2dxy;hessian(1,3)=df2dxz;
hessian(2,1)=df2dxy;hessian(2,2)=df2dy;hessian(2,3)=df2dyz;
hessian(3,1)=df2dxz;hessian(3,2)=df2dyz;hessian(3,3)=df2dz;

[xs1,xs2,xs3]=solve(dfdx,dfdy,dfdz,x,y,z);
xs1=vpa(xs1);xs2=vpa(xs2);xs3=vpa(xs3);

candidatos=double([xs1 xs2 xs3]);
fprintf('\n')
fprintf('CONDICAO DE 1A ORDEM\n')
confirm_grd=zeros(size(candidatos,1),3);
for i=1:size(candidatos,1)
confirm_grd(i,:)=double(round(subs(gradfx,[x,y,z],[candidatos(i,1),candidatos(i,2),candidatos(i,3)])));
fprintf('PONTO MINIMO RECORRENDO AS CONDICOES DE PRIMEIRA ORDEM:\n--> (%6.2f,%6.2f,%6.2f)\nCOM f=%6.3f\n ',candidatos(i,1),candidatos(i,2),candidatos(i,3),subs(thetas,[x,y,z],[candidatos(i,1),candidatos(i,2),candidatos(i,3)]))
end

fprintf('\n')
fprintf('CONDICAO DE 2A ORDEM\n')
confirm_hess=zeros(3,3,size(candidatos,1));k=0;
for i=1:size(candidatos,1)
confirm_hess(:,:,i)=double(round(subs(hessian,[x,y,z],[candidatos(i,1),candidatos(i,2),candidatos(i,3)])));
if confirm_hess(1,1,i)>=0 && confirm_hess(2,2,i)>=0 && confirm_hess(3,3,i)>=0 && det(confirm_hess(:,:,i))>=0
    fprintf('Hessiana de (%6.2f,%6.2f,%6.2f) e semi-positiva definida com determinante %6.3f\n',candidatos(i,1),candidatos(i,2),candidatos(i,3),det(confirm_hess(:,:,i)))
    k=k+1;minlocal(k,:)=candidatos(i,:);
else
    fprintf('Hessiana de (%6.2f,%6.2f,%6.2f) nao e semi-positiva definida\n',candidatos(i,1),candidatos(i,2),candidatos(i,3))
end
end

fprintf('\n')

fprintf('Minimos locais:\n')
for i=1:size(minlocal,1);fval(i)=subs(thetas,[x,y,z],[minlocal(i,1),minlocal(i,2),minlocal(i,3)]);fprintf('(%6.2f,%6.2f,%6.2f)-->f=%6.4f;\n',minlocal(i,1),minlocal(i,2),minlocal(i,3),fval(i));end

fprintf('\n')
confirm_hess=zeros(3,3,size(candidatos,1));k=0;
for i=1:size(candidatos,1)
confirm_hess(:,:,i)=double(round(subs(hessian,[x,y,z],[candidatos(i,1),candidatos(i,2),candidatos(i,3)])));
for j=1:3
    for jj=1:3
if isreal(confirm_hess(j,jj,i))
    k=k+1;
end
    end
end
end

if k==9
fprintf('Minimos globais:\n')
for i=1:size(minlocal,1);fval(i)=subs(thetas,[x,y,z],[minlocal(i,1),minlocal(i,2),minlocal(i,3)]);fprintf('(%6.2f,%6.2f,%6.2f)-->f=%6.4f;\n',minlocal(i,1),minlocal(i,2),minlocal(i,3),fval(i));end
end


end