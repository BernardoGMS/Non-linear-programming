function [x,y,D,lambda_,THETA,delta,ymais,ymenos]=prob2b


%% min((x1^2-x2)^2+2*(x2-x1)^4)

syms x1 x2 lambdaz
thetas=(x1^3-x2)^2+2*(x2-x1)^4; %funcao a minimizar

xp1=-.5:0.01:3.5;xp2=-.5:0.01:3.5;
[X1,X2]=meshgrid(xp1,xp2);

z=(X1.^3-X2).^2+2*(X2-X1).^4;

% Preparar o plot da funcao
titulo='Minimizacao - Metodo de Hooke-Jeeves';
            
scrsz = get(0,'ScreenSize');
posfig=[scrsz(3) scrsz(4) scrsz(3) scrsz(4)];

figure('Name',['Problema 2: ',titulo],...
'NumberTitle','off','OuterPosition',posfig);
hold on;

contour(X1,X2,z,[50:20:300 0:1:10],'ShowText','on');hold on;
xlabel('X1');ylabel('X2');

%% b) Resolver pelo Metodo de Hooke-Jeeves

%% METODO HOOKE-JEEVES

tdirect=1; %escolher linha de procura (1) ou passos discretos (2)

%% Passo de inicializa��o

epsilon=.01; %escalar maior que 0
x(1,:)=[0 3]; %ponto inicial de procura
d=[1 0;0 1];
y(1,:,1)=x;

if tdirect==1
    
    l=0.01;
    intervalo(1,1)=-10;intervalo(1,2)=10; %intervalo inicial de incerteza
    a(1)=intervalo(1,1);b(1)=intervalo(1,2);
    alpha1=0.618;
    lambda(1)=a(1)+(1-alpha1)*(b(1)-a(1));
    miu(1)=a(1)+alpha1*(b(1)-a(1));  
    
elseif tdirect==2
    
    delta(1)=.2;
    if delta(1)<epsilon
        error('delta tem de ser maior que epsilon')
    end
    ac_factor=2;
    ymais=[];ymenos=[];D=[];lambda_=[];
    
end

THETA(1)=subs(thetas,[x1,x2],[x(1,1),x(1,2)]);

k=1;
j=k;
n=2;
%% Passos principais

while j<n+1 % Passo 1

    if j==1;fprintf('Iteracao numero: %g, com x=(%6.3f,%6.3f)\n',k,x(end,1),x(end,2));end
    
    if tdirect==1

        lambdas=subs(thetas,[x1,x2],[y(j,1,k)+lambdaz*d(j,1),y(j,2,k)+lambdaz*d(j,2)]);
    lambdaj(j,k)=goldensection(lambdas,l,a,b,alpha1,lambda,miu);
    
    y(j+1,:,k)=y(j,:,k)+d(j,:)*lambdaj(j,k);
    
    plot(y(j:j+1,1,k),y(j:j+1,2,k),'ro-');hold on;
    
    if j<n
        j=j+1;
    else
    x(k+1,:)=y(n+1,:,k);
    if norm(x(k+1,:)-x(k,:))<epsilon 
        
     THETA(k+1)=subs(thetas,[x1,x2],[x(end,1),x(end,2)]);
     plot(x(:,1),x(:,2),'g*');hold on;plot(x(end,1),x(end,2),'go','LineWidth',5);hold on;
     for kl=1:size(x,1)
          ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
          set(ht,'HorizontalAlignment','right','FontSize',14)
     end
        ymais=0;ymenos=0;delta=0;
        if k==1;D=0;lambda_=0;end
     return
    else
     THETA(k+1)=subs(thetas,[x1,x2],[x(k+1,1),x(k+1,2)]);
     [y,j,k,DD,lambda__]=passo2(x,l,a,b,alpha1,lambda,miu,y,k,thetas); %Passo 2
     D(k-1,:)=DD;
     lambda_(k-1)=lambda__;
    end
    
    end

    elseif tdirect==2
        
        if (subs(thetas,[x1,x2],[y(j,1,k)+delta(k)*d(j,1),y(j,2,k)+delta(k)*d(j,2)]))<...
                (subs(thetas,[x1,x2],[y(j,1,k),y(j,2,k)]))
            
            ymais(j,:,k)=y(j,:,k)+delta(k)*d(j,:);
            fprintf('Sucesso! \n')
            
            y(j+1,:,k)=y(j,:,k)+d(j,:)*delta(k);
            
            plot(y(j:j+1,1,k),y(j:j+1,2,k),'ro-');hold on;

            [y,j,k,x,delta,THETA]=passo2_passos(x,y,k,ac_factor,j,delta,thetas,n,THETA,epsilon);
            
            if delta(k)<=epsilon
                
                plot(x(:,1),x(:,2),'g*');hold on;
                for kl=1:size(x,1)
                ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
                set(ht,'HorizontalAlignment','right','FontSize',14)
                end
                THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
                return
                
            end
            
        elseif (subs(thetas,[x1,x2],[y(j,1,k)+delta(k)*d(j,1),y(j,2,k)+delta(k)*d(j,2)]))>=...
                (subs(thetas,[x1,x2],[y(j,1,k),y(j,2,k)]))
            
            ymenos(j,:,k)=y(j,:,k)-delta(k)*d(j,:);
            fprintf('Fracasso! \n')
            
            y(j+1,:,k)=y(j,:,k)-d(j,:)*delta(k);
            
            plot(y(j:j+1,1,k),y(j:j+1,2,k),'ro-');hold on;
            
            [y,j,k,x,delta,THETA]=passo2_passos(x,y,k,ac_factor,j,delta,thetas,n,THETA,epsilon);
            
            if delta(k)<=epsilon
                
                plot(x(:,1),x(:,2),'g*');hold on;
                for kl=1:size(x,1)
                ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
                set(ht,'HorizontalAlignment','right','FontSize',14)
                end
                THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
                lambda_=0;D=0;
                return
                
            end
        end
    end
end

end

function [y,j,k,D,lambda_]=passo2(x,l,a,b,alpha1,lambda,miu,y,k,thetas)

    syms lambdaz x1 x2

        % Passo 2
    
    D=x(k+1,:)-x(k,:);
    
    lambdas=subs(thetas,[x1,x2],[x(k+1,1)+lambdaz*D(1),x(k+1,2)+lambdaz*D(2)]);
    
    lambda_=goldensection(lambdas,l,a,b,alpha1,lambda,miu);
    
    y(1,:,k+1)=x(k+1,:)+lambda_*D;
    
    yy(1,:)=y(3,:,k);yy(2,:)=y(1,:,k+1);plot(yy(:,1),yy(:,2),'r-.');hold on;
    
    j=1;k=k+1;

end

function [y,j,k,x,delta,THETA]=passo2_passos(x,y,k,ac_factor,j,delta,thetas,n,THETA,epsilon)

    syms x1 x2

        if j<n % Passo 2
                j=j+1;
            elseif j==n && (subs(thetas,[x1,x2],[y(n+1,1,k),y(n+1,2,k)]))<...
                (subs(thetas,[x1,x2],[x(k,1),x(k,2)])) % Passo 3
 
                [y,j,k,x,THETA,delta]=passo3_passos(x,y,k,ac_factor,n,thetas,delta,THETA);
            
            elseif j==n && (subs(thetas,[x1,x2],[y(n+1,1,k),y(n+1,2,k)]))>=...
                (subs(thetas,[x1,x2],[x(k,1),x(k,2)])) % Passo 4

            if delta<=epsilon

                THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
                return
                
            else
                
                delta(k+1)=delta(k)/2;
                
                y(1,:,k+1)=x(k,:);
                yy(1,:)=y(3,:,k);yy(2,:)=y(1,:,k+1);plot(yy(:,1),yy(:,2),'r-.');hold on
                x(k+1,:)=x(k,:);
                
                THETA(k+1)=subs(thetas,[x1,x2],[x(k+1,1),x(k+1,2)]);
                
                k=k+1;j=1;
                
            end
                
       end

end

function [y,j,k,x,THETA,delta]=passo3_passos(x,y,k,ac_factor,n,thetas,delta,THETA)

    syms x1 x2

        % Passo 3
    
    x(k+1,:)=y(n+1,:,k);
    
    THETA(k+1)=subs(thetas,[x1,x2],[x(k+1,1),x(k+1,2)]);

    y(1,:,k+1)=x(k+1,:)+ac_factor*(x(k+1,:)-x(k,:));
    
    yy(1,:)=y(3,:,k);yy(2,:)=y(1,:,k+1);plot(yy(:,1),yy(:,2),'r-.');hold on;
    
    delta(k+1)=delta(k);
    
    j=1;k=k+1;

end