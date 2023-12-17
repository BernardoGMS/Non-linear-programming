function [xx,yy,d,THETA,D,lambdaj]=prob5c


%% min((x1^2-x2)^2+2*(x2-x1)^4)

syms x y z lambdaz
thetas=prob5ab;
% thetas=(x1^3-x2)^2+2*(x2-x1)^4; %funcao a minimizar

gradiente=[diff(thetas,x),diff(thetas,y),diff(thetas,z)];

%% c) Resolver pelo Metodo de Davidon-Fletcher-Powell (DFP)

%% METODO DFP

%% Passo de inicialização

epsilon=0.01; %escalar maior que 0
xx(1,:)=[-10 -10 -10]; %ponto inicial de procura

yy(1,:,1)=xx;
D(:,:,1)=[1 0 0;0 1 0;0 0 1];

    l=0.001;
    intervalo(1,1)=0;intervalo(1,2)=5; %intervalo inicial de incerteza
    a(1)=intervalo(1,1);b(1)=intervalo(1,2); 
   
THETA(1)=subs(thetas,[x,y,z],[xx(1,1),xx(1,2),xx(1,3)]);

k=1;
j=k;
n=2;
%% Passos principais

while j<n+1 % Passo 1
    
    fprintf('-------------------------------- \n')

    if norm(subs(gradiente,[x,y,z],[yy(j,1,k),yy(j,2,k),yy(j,3,k)]))<epsilon
        
        THETA(k)=subs(thetas,[x,y,z],[xx(k,1),xx(k,2),xx(k,3)]);
          plot3(xx(:,1),xx(:,2),xx(:,3),'g*');hold on;
          for kl=1:size(xx,1)
          ht = text(xx(kl,1),xx(kl,2),xx(kl,2),['x_',int2str(kl)]);
          set(ht,'HorizontalAlignment','right')
          end
        return
        
    else

        d(j,:,k)=-D(:,:,k)*subs(gradiente,[x,y,z],[yy(j,1,k),yy(j,2,k),yy(j,3,k)])';

        lambdas=subs(thetas,[x,y,z],[yy(j,1,k)+lambdaz*d(j,1,k),yy(j,2,k)+lambdaz*d(j,2,k),yy(j,3,k)+lambdaz*d(j,3,k)]);
        lambdaj(j,k)=bisection(lambdas,l,a,b);
        
        fprintf('lambda %6.5f \n',lambdaj(j,k))
        
         yy(j+1,:,k)=yy(j,:,k)+d(j,:,k)*lambdaj(j,k);
         
         for jk=1:k
         plot3(yy(j:j+1,1,jk),yy(j:j+1,2,jk),yy(j:j+1,3,jk),'ro-');hold on;
         end
         
    if j<n % Passo 2
        
        p(j,:,k)=lambdaj(j,k)*d(j,:,k);

        q(j,:,k)=subs(gradiente,[x,y,z],[yy(j+1,1,k),yy(j+1,2,k),yy(j+1,3,k)])-subs(gradiente,[x,y,z],[yy(j,1,k),yy(j,2,k),yy(j,3,k)]);
        
        mult1=p(j,:,k)'*p(j,:,k);
        mult2=q(j,:,k)'*q(j,:,k);
        mult3=D(:,:,k)*q(j,:,k)'*q(j,:,k)*D(:,:,k);
        mult4=q(j,:,k)*D(:,:,k)*q(j,:,k)';
        D(:,:,k)=D(:,:,k)+mult1./mult2-mult3./mult4;
        
        grady=subs(gradiente,[x,y,z],[yy(j+1,1,k),yy(j+1,2,k),yy(j+1,3,k)]);grady2=subs(gradiente,[x,y,z],[yy(j,1,k),yy(j,2,k),yy(j,3,k)]);

        j=j+1;
        
    else
        
    yy(1,:,k+1)=yy(n+1,:,k);   
    xx(k+1,:)=yy(n+1,:,k);
    D(:,:,k+1)=[1 0 0;0 1 0;0 0 1];
    
    THETA(k)=subs(thetas,[x,y,z],[xx(k,1),xx(k,2),xx(k,3)]);
    k=k+1;
    j=1;
    
    end
        
    end

end

end
