function [x,y,d,THETA,alpha]=prob2d

%% min((x1^2-x2)^2+2*(x2-x1)^4)

syms x1 x2 lambdaz
thetas=(x1^3-x2)^2+2*(x2-x1)^4; %funcao a minimizar

gradiente=[diff(thetas,x1),diff(thetas,x2)];

xp1=-5:0.01:3.5;xp2=-5:0.01:3.5;
[X1,X2]=meshgrid(xp1,xp2);

z=(X1.^3-X2).^2+2*(X2-X1).^4;

% Preparar o plot da funcao
titulo='Minimizacao - Metodo dos GC de Fletcher and Reeves';
            
scrsz = get(0,'ScreenSize');
posfig=[scrsz(3) scrsz(4) scrsz(3) scrsz(4)];

figure('Name',['Problema 2: ',titulo],...
'NumberTitle','off','OuterPosition',posfig);
hold on;

contour(X1,X2,z,[100:500:4000 0:1:10],'ShowText','on');hold on;
xlabel('X1');ylabel('X2');

%% d) Resolver pelo Metodo dos Gradientes Conjugados de Fletcher and Reeves

%% METODO DE FLETCHER AND REEVES

%% Passo de inicialização

epsilon=0.05; %escalar maior que 0
x(1,:)=[-4 -5]; %ponto inicial de procura

y(1,:,1)=x;
d(:,:,1)=-[subs(gradiente(1),[x1,x2],[y(1,1,1),y(1,2,1)]) subs(gradiente(2),[x1,x2],[y(1,1,1),y(1,2,1)])];
    
    l=0.001;
    intervalo(1,1)=0;intervalo(1,2)=100; %intervalo inicial de incerteza
    a(1)=intervalo(1,1);b(1)=intervalo(1,2); 
   
THETA(1)=subs(thetas,[x1,x2],[x(1,1),x(1,2)]);

k=1;
j=k;
n=2;
%% Passos principais

while j<n+1 % Passo 1
    
    fprintf('-------------------------------- \n')
%     norm(subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)]))
    if norm(subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)]))<epsilon
        
        THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
          plot(x(:,1),x(:,2),'g*');hold on;
           for kl=1:size(x,1)
          ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
          set(ht,'HorizontalAlignment','right','FontSize',14)
           end
        return
        
    else
        
        lambdas=subs(thetas,[x1,x2],[y(j,1,k)+lambdaz*d(j,1,k),y(j,2,k)+lambdaz*d(j,2,k)]);
        lambdaj(j,k)=bisection(lambdas,l,a,b);
        
        fprintf('lambda %6.3f \n',lambdaj(j,k))
        
         y(j+1,:,k)=y(j,:,k)+d(j,:,k)*lambdaj(j,k);
         
         for jk=1:k
         plot(y(j:j+1,1,jk),y(j:j+1,2,jk),'ro-');hold on;
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
    
    THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
    k=k+1;
    j=1;
    
    end
        
    end

end

end
