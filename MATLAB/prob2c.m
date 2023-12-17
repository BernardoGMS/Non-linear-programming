function [x,y,d,lambda_,THETA,delta,ymais,ymenos]=prob2c

%% min((x1^2-x2)^2+2*(x2-x1)^4)

syms x1 x2 lambdaz
thetas=(x1^3-x2)^2+2*(x2-x1)^4; %funcao a minimizar

xp1=0:0.01:5;xp2=0:0.01:5;
[X1,X2]=meshgrid(xp1,xp2);

z=(X1.^3-X2).^2+2*(X2-X1).^4;

% Preparar o plot da funcao
titulo='Minimizacao - Metodo de Rosenbrock';
            
scrsz = get(0,'ScreenSize');
posfig=[scrsz(3) scrsz(4) scrsz(3) scrsz(4)];

figure('Name',['Problema 2: ',titulo],...
'NumberTitle','off','OuterPosition',posfig);
hold on;

contour(X1,X2,z,[0:2:20 100:1000:3000],'ShowText','on');hold on;
xlabel('X1');ylabel('X2');

%% c) Resolver pelo Metodo de Rosenbrock

%% METODO ROSENBROCK

tdirect=1; %escolher linha de procura (1) ou passos discretos (2)

%% Passo de inicialização

epsilon=0.05; %escalar maior que 0
x(1,:)=[3 2]; %ponto inicial de procura
d(:,:,1)=[1 0;0 1];
y(1,:,1)=x;

if tdirect==1
    
    l=0.1;
    intervalo(1,1)=-50;intervalo(1,2)=50; %intervalo inicial de incerteza
    a(1)=intervalo(1,1);b(1)=intervalo(1,2);
    alpha1=0.618;
    lambda(1)=a(1)+(1-alpha1)*(b(1)-a(1));
    miu(1)=a(1)+alpha1*(b(1)-a(1));  
    
elseif tdirect==2
    
    delta_inicial=.2;
    delta(1,1)=delta_inicial;delta(2,1)=delta(1,1);
    if delta(1)<0
        error('delta tem de ser positivo')
    end
    ex_factor=2; %expansion factor
    beta=-.7; %contraction factor
    
    if ex_factor<=1
        error('factor de expansao tem de ser maior que 1')
    end
    
    if beta<-1 || beta>0
        error('factor de contracao tem de estar no intervalo [-1;0]')
    end
    
end

THETA(1)=subs(thetas,[x1,x2],[x(1,1),x(1,2)]);

k=1;
j=k;
n=2;
%% Passos principais

while j<n+1 % Passo 1
    
    if tdirect==1
    fprintf('-------------------------------- \n')
    
    lambdas=subs(thetas,[x1,x2],[y(j,1,k)+lambdaz*d(j,1,k),y(j,2,k)+lambdaz*d(j,2,k)]);
    lambdaj(j,k)=goldensection(lambdas,l,a,b,alpha1,lambda,miu);
    
    y(j+1,:,k)=y(j,:,k)+d(j,:,k)*lambdaj(j,k);
    
    plot(y(j:j+1,1,k),y(j:j+1,2,k),'ro-');hold on;
    
    if j<n
        j=j+1;
    else
    x(k+1,:)=y(n+1,:,k);
%     norm(x(k+1,:)-x(k,:))
    if norm(x(k+1,:)-x(k,:))<epsilon 
        
     THETA(k+1)=subs(thetas,[x1,x2],[x(end,1),x(end,2)]);
     plot(x(:,1),x(:,2),'g*');hold on;
     for kl=1:size(x,1)
          ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
          set(ht,'HorizontalAlignment','right','FontSize',14)
     end
     ymais=0;ymenos=0;delta=0;lambda_=0;
     return
    else
     THETA(k+1)=subs(thetas,[x1,x2],[x(k+1,1),x(k+1,2)]);
     [y,j,k,d]=passo2_rosenbrock(x,y,k,lambdaj,d,n); %Passo 2

    end
    
    end
    
    elseif tdirect==2

        if (subs(thetas,[x1,x2],[y(j,1,k)+delta(j,k)*d(j,1,k),y(j,2,k)+delta(j,k)*d(j,2,k)]))<...
                (subs(thetas,[x1,x2],[y(j,1,k),y(j,2,k)]))
            
            fprintf('Sucesso com j=%g e k=%g! \n',j,k)
            
            y(j+1,:,k)=y(j,:,k)+d(j,:,k)*delta(j,k);
            delta(j,k)=delta(j,k)*ex_factor;
            fprintf('delta(%6.2f) \n',delta(j,k))
            
            fprintf('y=(%6.2f;%6.2f) \n',y(j,1,k),y(j,2,k))
            fprintf('y+deltaj*dj=(%6.2f;%6.2f) \n',y(j+1,1,k),y(j+1,2,k))
            fprintf('d=(%6.2f;%6.2f) \n',d(j,1,k),d(j,2,k))
            
            plot(y(j:j+1,1,k),y(j:j+1,2,k),'ro-');hold on;

            ver_k=k;indicador=0;

            [y,j,k,x,delta,THETA,d,indicador]=passo1_passos(x,y,k,j,delta,thetas,n,THETA,epsilon,d,delta_inicial,indicador);
            
            if ver_k==k 
                
            if size(x,1)>ver_k

                if indicador==1
                plot(x(:,1),x(:,2),'b*');hold on;
                for kl=1:size(x,1)
                ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
                set(ht,'HorizontalAlignment','right','FontSize',14)
                end
                for jk=1:k
                    plot(y(:,1,jk),y(:,2,jk),'r.');hold on;
                end
                lambda_=0;ymais=0;ymenos=0;
                return
                
                end
            
            end

                if indicador==1
                lambda_=0;ymais=0;ymenos=0;
                THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
                plot(x(:,1),x(:,2),'b*');hold on;
                for kl=1:size(x,1)
                ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
                set(ht,'HorizontalAlignment','right','FontSize',14)
                end
                for jk=1:k
                    plot(y(:,1,jk),y(:,2,jk),'r.');hold on;
                end
                return
                end

          
            end

        elseif (subs(thetas,[x1,x2],[y(j,1,k)+delta(j,k)*d(j,1,k),y(j,2,k)+delta(j,k)*d(j,2,k)]))>=...
                (subs(thetas,[x1,x2],[y(j,1,k),y(j,2,k)]))
            
            fprintf('Fracasso com j=%g e k=%g! \n',j,k)
            
            fprintf('y=(%6.2f;%6.2f) \n',y(j,1,k),y(j,2,k))
            fprintf('y+delta*dj =(%6.2f;%6.2f)\n',y(j,1,k)+d(j,1,k)*delta(j,k),y(j,2,k)+d(j,2,k)*delta(j,k))
            fprintf('d=(%6.2f;%6.2f) \n',d(j,1,k),d(j,2,k))
            
            y(j+1,:,k)=y(j,:,k);
            delta(j,k)=beta*delta(j,k);
            fprintf('delta(%6.2f) \n',delta(j,k))

            plot(y(j:j+1,1,k),y(j:j+1,2,k),'ro-');hold on;
            ver_k=k;indicador=0;

            [y,j,k,x,delta,THETA,d,indicador]=passo1_passos(x,y,k,j,delta,thetas,n,THETA,epsilon,d,delta_inicial,indicador);
            
            if ver_k==k 

                if indicador==1
                THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
                plot(x(:,1),x(:,2),'b*');hold on;
                for kl=1:size(x,1)
                ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
                set(ht,'HorizontalAlignment','right','FontSize',14)
                end
                for jk=1:k
                    plot(y(:,1,jk),y(:,2,jk),'r.');hold on;
                end
                return
                
                end
            
                if indicador==1 
                fprintf('abs(delta(1,k))<=epsilon && abs(delta(2,k))<=epsilon\n')

                THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
                plot(x(:,1),x(:,2),'b*');hold on;
                for kl=1:size(x,1)
                ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
                set(ht,'HorizontalAlignment','right','FontSize',14)
                end
                for jk=1:k
                    plot(y(:,1,jk),y(:,2,jk),'r.');hold on;
                end
                return
                end
          
            end

        end
        
        lambda_=0;ymais=0;ymenos=0;
        
        if isnan((subs(thetas,[x1,x2],[y(j,1,k)+delta(j,k)*d(j,1,k),y(j,2,k)+delta(j,k)*d(j,2,k)]))) && k>1
        fprintf('Processo viciado');
        return
        end
    
    end
    
    
end

end

function [y,j,k,d]=passo2_rosenbrock(x,y,k,lambdaj,d,n)

    fprintf('Passo 2 com k= %g\n',k)

        % Passo 2
    
    y(1,:,k+1)=x(k+1,:);
    
    yy(1,:)=y(3,:,k);yy(2,:)=y(1,:,k+1);plot(yy(:,1),yy(:,2),'r-.');hold on;
    
    j=1;k=k+1;
    
        %Passo 3 - Processo de Gram-Schmit

    a=zeros(n,n);

   for kk=1:n
    for jj=kk:n
    a(kk,:)=a(kk,:)+lambdaj(jj,k-1)*d(jj,:,k-1);
    end
   end

    b(1,:)=a(1,:);
    d_(1,:)=b(1,:)/norm(b(1,:));

    for jj=2:n
        for ii=1:jj-1
         b(jj,:)=a(jj,:)-a(jj,:)*d_(jj-1,:)'*d_(jj-1,:);
        end
        d_(jj,:)=b(jj,:)/norm(b(jj,:));
    end
    
    d(:,:,k)=d_;
    

end

function [y,j,k,x,delta,THETA,d,indicador]=passo1_passos(x,y,k,j,delta,thetas,n,THETA,epsilon,d,delta_inicial,indicador)

    syms x1 x2

    fprintf('Passo 1 com k= %g\n',k)
        
   
        if j<n % Passo 2
                j=j+1;
        elseif j==n && (subs(thetas,[x1,x2],[y(n+1,1,k),y(n+1,2,k)]))<...
                (subs(thetas,[x1,x2],[y(1,1,k),y(1,2,k)])) % se cada um dos trials do passo 1 foram sucesso
                    
                    fprintf('se cada um dos trials do passo 1 foram sucesso\n')
                    y(1,:,k)=y(n+1,:,k);
                    j=1;
            
        elseif j==n && (subs(thetas,[x1,x2],[y(n+1,1,k),y(n+1,2,k)]))==...
                (subs(thetas,[x1,x2],[y(1,1,k),y(1,2,k)])) % quando cada um dos ultimos trials do passo 1 foram falha
            
            fprintf('quando cada um dos ultimos trials do passo 1 foram falha\n')
            
            if (subs(thetas,[x1,x2],[y(n+1,1,k),y(n+1,2,k)]))<...
                (subs(thetas,[x1,x2],[x(k,1),x(k,2)])) % se pelo menos um sucesso encontrado no passo 1 durante iteracao k
            
            fprintf('se pelo menos um sucesso encontrado no passo 1 durante iteracao k\n')
            
            %passo 3
            [y,j,k,x,THETA,delta,d,indicador]=passo3_passos(x,y,k,n,thetas,delta,THETA,d,epsilon,delta_inicial,indicador);
            
            if indicador==1
     
                fprintf('abs(x(k,:)-x(k-1,:))<epsilon\n')
                plot(x(:,1),x(:,2),'b*');hold on;
                for kl=1:size(x,1)
                ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
                set(ht,'HorizontalAlignment','right','FontSize',14)
                end
                for jk=1:k
                    plot(y(:,1,jk),y(:,2,jk),'r.');hold on;
                end
                return
   
            end
            
            elseif j==n && (subs(thetas,[x1,x2],[y(n+1,1,k),y(n+1,2,k)]))==...
                (subs(thetas,[x1,x2],[x(k,1),x(k,2)])) % nenhum sucesso encontrado no passo 1 durante iteracao k
                
            fprintf('nenhum sucesso encontrado no passo 1 durante iteracao k\n')
                %abs(delta(1,k))<=epsilon &&
            if  abs(delta(2,k))<=epsilon
                plot(x(:,1),x(:,2),'b*');hold on;
                for kl=1:size(x,1)
                ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
                set(ht,'HorizontalAlignment','right','FontSize',14)
                end
                for jk=1:k
                    plot(y(:,1,jk),y(:,2,jk),'r.');hold on;
                end
                indicador=1;
                fprintf('abs(delta(1,k))<=epsilon && abs(delta(2,k))<=epsilon\n')
                THETA(k)=subs(thetas,[x1,x2],[x(k,1),x(k,2)]);
                return
                
            else
                
                y(1,:,k)=y(n+1,:,k);
                
                j=1;
                
            end

            end
       end
end

function [y,j,k,x,THETA,delta,d,indicador]=passo3_passos(x,y,k,n,thetas,delta,THETA,d,epsilon,delta_inicial,indicador)

    syms x1 x2
       
        fprintf('Passo 3 com k= %g\n',k)

        % Passo 3
    
    x(k+1,:)=y(n+1,:,k);
    
    if abs(norm(x(k+1,:)-x(k,:)))<epsilon
        
        indicador=1;
        fprintf('abs(x(k+1,:)-x(k,:))<epsilon\n')
        plot(x(:,1),x(:,2),'b*');hold on;
        for jk=1:k
           plot(y(:,1,jk),y(:,2,jk),'r.');hold on;
        end
        THETA(k+1)=subs(thetas,[x1,x2],[x(k+1,1),x(k+1,2)]);j=n;
                return
                
    else

    syms lambda1 lambda2

    lamb1=[lambda1,lambda2];

    deltax=x(k+1,:)-x(k,:);

    equ=0;
    
    for i=1:n
   
        equ=equ+lamb1(i)*d(i,:,k);
    
    end

    [d1,d2]=solve(deltax(1)-equ(1),deltax(2)-equ(2));
    
    %Passo 3 - Processo de Gram-Schmit

    lambdaj(:,k)=[d1; d2]; 
        
    a=zeros(n,n);

   for kk=1:n
    for jj=kk:n
    a(kk,:)=a(kk,:)+lambdaj(jj,k)*d(jj,:,k);
    end
   end

    b(1,:)=a(1,:);
    d_(1,:)=b(1,:)/norm(b(1,:));

    for jj=2:n
        for ii=1:jj-1
         b(jj,:)=a(jj,:)-a(jj,:)*d_(jj-1,:)'*d_(jj-1,:);
        end
        d_(jj,:)=b(jj,:)/norm(b(jj,:));
    end
    
    d(:,:,k+1)=d_;
    
    delta(1,k+1)=delta_inicial;delta(2,k+1)=delta(1,k+1);
    
    THETA(k+1)=subs(thetas,[x1,x2],[x(k+1,1),x(k+1,2)]);
 
    y(1,:,k+1)=x(k+1,:);
   
    yy(1,:)=y(3,:,k);yy(2,:)=y(1,:,k+1);plot(yy(:,1),yy(:,2),'r-.');hold on;

    j=1;k=k+1;
        
    end
     

end
