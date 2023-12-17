function prob2a

%% min((x1^2-x2)^2+2*(x2-x1)^4)

syms x1 x2 real
syms lambdaz
thetas=(x1^3-x2)^2+2*(x2-x1)^4; %funcao a minimizar
%% a) Condicoes necessarias de minimo:

% Condicao de 1ordem para x* ser um optimo local:
% E d tal que grad(f(x*))t.d<0 --> E delta>0 tal que f(x*+lambda.d)<f(x*)
%para cada lambda pertencente [0, delta];

gradfx=sym(zeros(2,1));hessian=sym(zeros(2,2));
dfdx1=diff(thetas,x1);dfdx2=diff(thetas,x2);
df2dx1=diff(dfdx1,x1);df2dx2=diff(dfdx2,x2);df2dx1x2=diff(dfdx2,x1);
gradfx(1,1)=dfdx1;gradfx(2,1)=dfdx2;
hessian(1,1)=df2dx1;hessian(2,1)=df2dx1x2;hessian(1,2)=df2dx1x2;hessian(2,2)=df2dx2;

[xs1,xs2]=solve(dfdx1,dfdx2,x1,x2);
xs1=vpa(xs1);xs2=vpa(xs2);

candidatos=double([xs1 xs2]);

% Condicao de 1a ordem
confirm_grd=zeros(size(candidatos,1),2);
for i=1:size(candidatos,1)
confirm_grd(i,:)=double(round(subs(gradfx,[x1,x2],[candidatos(i,1),candidatos(i,2)])));
fprintf('%g --> (%6.2f,%6.2f)\n',i,candidatos(i,1),candidatos(i,2))
end

fprintf('\n')
if confirm_grd==zeros(size(confirm_grd));fprintf('Todos os candidatos a minimo local satisfazem a condicao de minimo de 1a ordem!\n');end
fprintf('\n')

% Condicao de 2a ordem
confirm_hess=zeros(2,2,size(candidatos,1));k=0;
for i=1:size(candidatos,1)
confirm_hess(:,:,i)=double(round(subs(hessian,[x1,x2],[candidatos(i,1),candidatos(i,2)])));
if confirm_hess(1,1,i)>=0 && confirm_hess(2,2,i)>=0 && confirm_hess(1,1,i)*confirm_hess(2,2,i)-confirm_hess(1,2,i)^2>=0
    fprintf('Hessiana de (%6.2f,%6.2f) e semi-positiva definida\n',candidatos(i,1),candidatos(i,2))
    k=k+1;minlocal(k,:)=candidatos(i,:);
else
    fprintf('Hessiana de (%6.2f,%6.2f) nao e semi-positiva definida\n',candidatos(i,1),candidatos(i,2))
end
end

fprintf('\n')

fprintf('Minimos locais:\n')
for i=1:size(minlocal,1);fval(i)=subs(thetas,[x1,x2],[minlocal(i,1),minlocal(i,2)]);fprintf('(%6.2f,%6.2f)-->f=%6.4f;\n',minlocal(i,1),minlocal(i,2),fval(i));end

xp1=-3:0.01:3;xp2=-3:0.01:3;
[X1,X2]=meshgrid(xp1,xp2);

z=(X1.^3-X2).^2+2*(X2-X1).^4;

% Preparar o plot da funcao
titulo='Isolinhas da funcao a minimizar';
            
scrsz = get(0,'ScreenSize');
posfig=[scrsz(3) scrsz(4) scrsz(3) scrsz(4)];

figure('Name',['Problema 2: ',titulo],...
'NumberTitle','off','OuterPosition',posfig);
hold on;

contour(X1,X2,z,[0:.4:5],'ShowText','on');hold on;
xlabel('X1');ylabel('X2');grid on

% Condicao de 2ordem para x* ser um optimo local:
% f:En-->E1 é duplamente diff em x*, se x* e um minimo local, entao
% grad(f(x*))=0 e H(x*) e positivo semi-definido


end
