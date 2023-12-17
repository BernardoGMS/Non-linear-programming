function [x,y,d,THETA,D,lambdaj]=prob2e


%% min((x1^2-x2)^2+2*(x2-x1)^4)

syms x1 x2 lambdaz
thetas=(x1^3-x2)^2+2*(x2-x1)^4; %funcao a minimizar

gradiente=[diff(thetas,x1),diff(thetas,x2)];

xp1=-5:0.01:3.5;xp2=-5:0.01:3.5;
[X1,X2]=meshgrid(xp1,xp2);

z=(X1.^3-X2).^2+2*(X2-X1).^4;

% Preparar o plot da funcao
titulo='Minimizacao - Metodo de Davidon-Fletcher-Powell';
            
scrsz = get(0,'ScreenSize');
posfig=[scrsz(3) scrsz(4) scrsz(3) scrsz(4)];

figure('Name',['Problema 2: ',titulo],...
'NumberTitle','off','OuterPosition',posfig);
hold on;

contour(X1,X2,z,[100:500:4000 0:1:10],'ShowText','on');hold on;
xlabel('X1');ylabel('X2');

%% e) Resolver pelo Metodo de Davidon-Fletcher-Powell (DFP)

%% METODO DFP

%% Passo de inicialização

epsilon=0.05; %escalar maior que 0
x(1,:)=[-4 -5]; %ponto inicial de procura

y(1,:,1)=x;
D(:,:,1)=[1 0;0 1];

    l=0.001;
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
          plot(x(:,1),x(:,2),'g*');hold on;
          for kl=1:size(x,1)
          ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
          set(ht,'HorizontalAlignment','right','FontSize',14)
          end
        return
        
    else

        d(j,:,k)=-D(:,:,k)*subs(gradiente,[x1,x2],[y(j,1,k),y(j,2,k)])';

        lambdas=subs(thetas,[x1,x2],[y(j,1,k)+lambdaz*d(j,1,k),y(j,2,k)+lambdaz*d(j,2,k)]);
        lambdaj(j,k)=bisection(lambdas,l,a,b);
        
        fprintf('lambda %6.5f \n',lambdaj(j,k))
        
         y(j+1,:,k)=y(j,:,k)+d(j,:,k)*lambdaj(j,k);
         
         for jk=1:k
         plot(y(j:j+1,1,jk),y(j:j+1,2,jk),'ro-');hold on;
         end
         
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
