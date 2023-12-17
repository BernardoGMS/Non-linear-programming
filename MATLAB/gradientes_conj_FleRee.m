function [x,TH]=gradientes_conj_FleRee(thetas)

%% min((x1^2-x2)^2+2*(x2-x1)^4)

syms x1 x2 lambdaz

gradiente=[diff(thetas,x1),diff(thetas,x2)];
[xs1,xs2]=solve(gradiente(1),gradiente(2));
xs1=double(xs1);xs2=double(xs2);ii=0;
for i=1:size(xs1);if isreal(xs1(i)) && isreal(xs2(i));ii=ii+1;xxs1(ii,1)=xs1(i);xxs2(ii,1)=xs2(i);fval(ii,1)=subs(thetas,[x1,x2],[xxs1(ii,1),xxs2(ii,1)]);end; end
indd=find(fval==min(min(fval)));
x(1,:)=[xxs1(indd,1),xxs2(indd,1)];

%% Resolver pelo Metodo dos Gradientes Conjugados de Fletcher and Reeves

%% METODO DE FLETCHER AND REEVES

%% Passo de inicialização

epsilon=0.05; %escalar maior que 0

y(1,:,1)=x(1,:);
d(:,:,1)=-[subs(gradiente(1),[x1,x2],[y(1,1,1),y(1,2,1)]) subs(gradiente(2),[x1,x2],[y(1,1,1),y(1,2,1)])];
    
    l=0.01;
    intervalo(1,1)=-20;intervalo(1,2)=20; %intervalo inicial de incerteza
    a(1)=intervalo(1,1);b(1)=intervalo(1,2); 
   
THETA(1)=subs(thetas,[x1,x2],[x(1,1),x(1,2)]);

k=1;
j=k;
n=2;
%% Passos principais

while j<n+1 % Passo 1
    
    fprintf('-------------------------------- \n')
    
    if norm(subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)]))<epsilon
        
        THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);TH=THETA(end);
          plot(x(:,1),x(:,2),'ro-');hold on;
     
        return
        
    else
        
        lambdas=subs(thetas,[x1,x2],[y(j,1,k)+lambdaz*d(j,1,k),y(j,2,k)+lambdaz*d(j,2,k)]);
        lambdaj(j,k)=bisection(lambdas,l,a,b);
        
        fprintf('lambda %6.3f \n',lambdaj(j,k))
        
         y(j+1,:,k)=y(j,:,k)+d(j,:,k)*lambdaj(j,k);
         
         for jk=1:k
         plot(y(j:j+1,1,jk),y(j:j+1,2,jk),'g*-');hold on;
         end
         
    if j<n % Passo 2
        
        alpha(j,k)=norm(subs(gradiente,[x1,x2],[y(j+1,1,k),y(j+1,2,k)]))^2/norm(subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)]))^2;
        d(j+1,:,k)=-subs(gradiente,[x1,x2],[y(j+1,1,k),y(j+1,2,k)])+alpha(j,k)*d(j,:,k);
        
        fprintf('alpha %6.3f \n',alpha(j,k))
        grady=subs(gradiente,[x1,x2],[y(j+1,1,k),y(j+1,2,k)]);grady2=subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)]);
        fprintf('grad yj (%6.3f,%6.3f) \n',grady2(1,1),grady2(1,2))
        fprintf('grad yj+1 (%6.3f,%6.3f) \n',grady(1,1),grady(1,2))
        fprintf('dj (%6.3f,%6.3f) \n',d(j,1,k),d(j,2,k))
        fprintf('dj+1 (%6.3f,%6.3f) \n',d(j+1,1,k),d(j+1,2,k))

        j=j+1;
    else % Passo 3
        
    y(1,:,k+1)=y(n+1,:,k);   
    x(k+1,:)=y(n+1,:,k);
    d(1,:,k+1)=-subs(gradiente,[x1,x2],[y(1,1,k+1),y(1,2,k+1)]);
    
    THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);TH=THETA(end);
    k=k+1;
    j=1;
    
    end
        
    end

end



end