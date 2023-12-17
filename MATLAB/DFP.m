function [x,THETA]=DFP(thetas)

%% min((x1^2-x2)^2+2*(x2-x1)^4)

syms x1 x2 lambdaz

gradiente=[diff(thetas,x1),diff(thetas,x2)];
[xs1,xs2]=solve(gradiente(1),gradiente(2));
xs1=double(xs1);xs2=double(xs2);ii=0;
for i=1:size(xs1);if isreal(xs1(i)) && isreal(xs2(i));ii=ii+1;xxs1(ii,1)=xs1(i);xxs2(ii,1)=xs2(i);fval(ii,1)=subs(thetas,[x1,x2],[xxs1(ii,1),xxs2(ii,1)]);end; end
indd=find(fval==min(min(fval)));
x(1,:)=[xxs1(indd,1),xxs2(indd,1)];

%% Resolver pelo Metodo de Davidon-Fletcher-Powell (DFP)

%% METODO DFP

%% Passo de inicialização

epsilon=0.01; %escalar maior que 0

y(1,:,1)=x;
D(:,:,1)=[1 0;0 1];

    l=0.01;
    intervalo(1,1)=0;intervalo(1,2)=5; %intervalo inicial de incerteza
    a(1)=intervalo(1,1);b(1)=intervalo(1,2); 
   
THETA(1)=subs(thetas,[x1,x2],[x(1,1),x(1,2)]);

k=1;
j=k;
n=2;
%% Passos principais

while j<n+1 % Passo 1
    
    fprintf('-------------------------------- \n')
    
    if norm(subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)]))<epsilon
        
        THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
        return
        
    else

        d(j,:,k)=-D(:,:,k)*subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)])';

        lambdas=subs(thetas,[x1,x2],[y(j,1,k)+lambdaz*d(j,1,k),y(j,2,k)+lambdaz*d(j,2,k)]);
        lambdaj(j,k)=bisection(lambdas,l,a,b);
        
        fprintf('lambda %6.5f \n',lambdaj(j,k))
        
         y(j+1,:,k)=y(j,:,k)+d(j,:,k)*lambdaj(j,k);
         
    if j<n % Passo 2
        
        p(j,:,k)=lambdaj(j,k)*d(j,:,k);

        q(j,:,k)=subs(gradiente,[x1,x2],[y(j+1,1,k),y(j+1,2,k)])-subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)]);

        mult1=p(j,:,k)'*p(j,:,k);
        mult2=q(j,:,k)'*q(j,:,k);
        mult3=D(:,:,k)*q(j,:,k)'*q(j,:,k)*D(:,:,k);
        mult4=q(j,:,k)*D(:,:,k)*q(j,:,k)';
        D(:,:,k)=D(:,:,k)+mult1./mult2-mult3./mult4;
        
        grady=subs(gradiente,[x1,x2],[y(j+1,1,k),y(j+1,2,k)]);grady2=subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)]);

        j=j+1;
        
    else
        
    y(1,:,k+1)=y(n+1,:,k);   
    x(k+1,:)=y(n+1,:,k);
    D(:,:,k+1)=[1 0 ;0 1];
    
    THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
    k=k+1;
    j=1;
    
    end
        
    end

end

end