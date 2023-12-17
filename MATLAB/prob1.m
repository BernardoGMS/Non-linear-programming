function [intervalo,lambda,miu,theta,n,compincert]=prob1

%% min(e^-lambda+lambda^2)

syms lambdas
fprintf('------------------- PROBLEMA 1 -------------------\n')
fprintf('              Funcao a minimizar\n')
thetas=exp(-lambdas)+lambdas^2; %funcao a minimizar
fprintf('           ');disp(thetas)
fprintf('--------------------------------------------------\n')
%% Verificar os zeros pela derivada de primeira ordem da funcao

diffthetas=diff(thetas,lambdas);
lambda_zero=double(solve(diffthetas,lambdas));

otimo=subs(thetas,lambdas,lambda_zero);
fprintf('Zero da funcao em lambda (4 digitos) = %6.4f\n',lambda_zero)
fprintf('Minimo da funcao: (4 digitos) = %6.7f\n',otimo)

%% Passo de inicializa��o

l=0.1; %comprimento final permitido de incerteza (l>0)
intervalo(1,1)=0;intervalo(1,2)=5; %intervalo inicial de incerteza
compincert(1)=abs(intervalo(1,1)-intervalo(1,2));
a(1)=intervalo(1,1);b(1)=intervalo(1,2);
alpha1=0.618;
lambda(1)=a(1)+(1-alpha1)*(b(1)-a(1));
miu(1)=a(1)+alpha1*(b(1)-a(1));

theta(1,1)=subs(thetas,lambdas,lambda(1));
theta(1,2)=subs(thetas,lambdas,miu(1));

% Preparar o plot da funcao
titulo='Minimizacao - Metodo da Seccao Dourada';
            
scrsz = get(0,'ScreenSize');
posfig=[scrsz(3) scrsz(4) scrsz(3) scrsz(4)];

figure('Name',['Problema 1: ',titulo],...
'NumberTitle','off','OuterPosition',posfig);
hold on;

X1=a(1):.01:b(1);THETAPLOT=exp(-X1)+X1.^2;DTHETAPLOT=-exp(-X1)+2*X1;
minim=solve(-exp(-lambdas)+2*lambdas,lambdas);
plot(X1,THETAPLOT,'b-');hold on; % Plot da funcao no intervalo inicial
plot(X1,DTHETAPLOT,'g-');hold on; % Plot da derivada funcao no intervalo inicial
line([minim minim],[-2 15],'Color','g','LineWidth',2,'LineStyle','--');line([0 5],[0 0],'Color','k','LineWidth',2,'LineStyle','-');grid on
plot(lambda(1),theta(1),'ro');hold on; % Plot das primeiras solucoes
plot(miu(1),theta(2),'ro');hold on;plot(minim,0,'Color','g','LineWidth',4,'LineStyle','o');
xlabel('\lambda');ylabel('\theta(\lambda)');
h=legend('\theta(\lambda)','d\theta(\lambda)/d\lambda','d\theta(\lambda)/d\lambda=0','Location','NorthEastOutside');set(h);

%% Passos principais

k=1;
while b(k)-a(k)>l
    
    if theta(k,1)>theta(k,2) % Passo 2
        
        a(k+1)=lambda(k);
        b(k+1)=b(k);
        lambda(k+1)=miu(k);
        miu(k+1)=a(k+1)+alpha1*(b(k+1)-a(k+1));
        
        theta(k+1,1)=theta(k,2);
        theta(k+1,2)=subs(thetas,lambdas,miu(k+1));
        
        plot(miu(k:k+1),theta(k:k+1,2),'r>-');hold on
        ht = text(miu(k+1),theta(k+1,2),['  ',int2str(k),'a']);
        set(ht,'HorizontalAlignment','right','FontSize',14)
        
    else % Passo 3
        
        a(k+1)=a(k);
        b(k+1)=miu(k);
        miu(k+1)=lambda(k);
        lambda(k+1)=a(k+1)+(1-alpha1)*(b(k+1)-a(k+1));
        
        theta(k+1,1)=subs(thetas,lambdas,lambda(k+1));
        theta(k+1,2)=theta(k,1);
        
         plot(lambda(k:k+1),theta(k:k+1,1),'r<-');hold on
         ht = text(miu(k+1),theta(k+1,2),['  ',int2str(k),'a']);
         set(ht,'HorizontalAlignment','right','FontSize',14)
        
    end
    
    k=k+1; % Passo 4
    compincert(k)=abs(b(k)-a(k)); %comprimento do intervalo de incerteza
    intervalo(k,:)=[a(k) b(k)]; %intervalo de incerteza

end

n=ceil(log((b(end)-a(end))/(b(1)-a(1)))/log(alpha1)+1);
theta(end,:)=[];lambda(end)=[];miu(end)=[];

end