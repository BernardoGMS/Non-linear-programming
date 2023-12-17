function [x,VIOLG,VIOLI,miu,miv,theta,THETA]=prob4f(metodo)

syms x1 x2

% Inicializacao
[x,miu,miv,thetaz,igthetaz,desigthetaz,alpha]=prob4d;
Nig=size(igthetaz,1);Ndes=size(desigthetaz,1);

theta(1,1)=subs(thetaz,[x1,x2],[x(1,1),x(1,2)]);

epsilon=10e-5;
VIOLG=zeros(1,Ndes);VIOLI=zeros(1,Nig);

k=1;
VIOLG(k,:)=max(subs(desigthetaz,[x1,x2],[x(k,1),x(k,2)]),0);
VIOLI(k,:)=max(abs(subs(igthetaz,[x1,x2],[x(k,1),x(k,2)])));

% Primeira vez:
% Rotina interna (Minimizacao da funcao penalizada - sem constrangimentos)
if strcmp(metodo,'GC')
[xx,THETA(k,1)]=gradientes_conj_FleRee(alpha);
elseif strcmp(metodo,'RS')
[xx,THETA(k,1)]=Rosenbrock(alpha);
else
error('Metodo desconhecido!!!!')
end

for i=1:25

x(k+1,:)=xx(end,:);
theta(k+1,1)=subs(thetaz,[x1,x2],[x(k+1,1),x(k+1,2)]);

plot(x(k:k+1,1),x(k:k+1,2),'g*-');hold on;


VIOLG(k+1,:)=max(subs(desigthetaz,[x1,x2],[x(k+1,1),x(k+1,2)]),0);
VIOLI(k+1,:)=max(abs(subs(igthetaz,[x1,x2],[x(k+1,1),x(k+1,2)])));
activeG(1,:)=zeros(1,Ndes);activeI(1,:)=zeros(1,Ndes);

%% Verificacao das violacoes das desigualdades
for j=1:Ndes
if VIOLG(k+1,j)>epsilon && activeG(k,j)==0
    activeG(k+1,j)=1;
    miu(k+1,j)=.1;fprintf('Activa %6.2f miu(%s,%s)...    ',miu(k+1,j),int2str(k+1),int2str(j))
elseif VIOLG(k+1,j)>epsilon && activeG(k,j)==1
    activeG(k+1,j)=1;
    miu(k+1,j)=10*miu(k,j);fprintf('Actualiza %6.2f miu(%s,%s)...    ',miu(k+1,j),int2str(k+1),int2str(j)) 
else activeG(k+1,j)=0;miu(k+1,j)=miu(k,j);fprintf('Nao actualiza %6.2f miu(%s,%s)...    ',miu(k+1,j),int2str(k+1),int2str(j))
end
fprintf('\n')
end

%% Verificacao das violacoes das igualdades

for j=1:Nig
if VIOLI(k+1,j)>epsilon && activeI(k,j)==0
    activeI(k+1,j)=1;
    miv(k+1,j)=.1;fprintf('Activa %6.2f miu(%s,%s)...    ',miv(k+1,j),int2str(k+1),int2str(j))
elseif VIOLI(k+1,j)>epsilon&& activeI(k,j)==1
    activeI(k+1,j)=1;
    miv(k+1,j)=10*miv(k,j);fprintf('Actualiza %6.2f miu(%s,%s)...    ',miv(k+1,j),int2str(k+1),int2str(j)) 
else activeI(k+1,j)=0;miv(k+1,j)=miv(k,j);fprintf('Nao actualiza %6.2f miu(%s,%s)...    ',miv(k+1,j),int2str(k+1),int2str(j))
end
fprintf('\n')
end

    k=k+1;
    
    alpha=thetaz+miu(k,:)*desigthetaz.^2+miv(k,:)*igthetaz.^2;
    fprintf('\nFUNCAO DESCONSTRANGIDA F(xk,miu1,miu2):\n')
    pretty(vpa(alpha))
    
    kk=0;
    for j=1:Ndes
        if VIOLG(k,j)<epsilon
        kk=kk+1;
        end
    end
    
    for j=1:Nig
        if VIOLI(k,j)<epsilon
        kk=kk+1;
        end
    end
    
    pretty(alpha)
    
%% Procurar nova solucao
    if strcmp(metodo,'GC')
    [xx,THETA(k,1)]=gradientes_conj_FleRee(alpha);
    elseif strcmp(metodo,'RS')
    [xx,THETA(k,1)]=Rosenbrock(alpha);
    else
    error('Metodo desconhecido!!!!')
    end

    if k>2
    if abs(THETA(k,:)-THETA(k-1,:))<epsilon && kk==Nig+Ndes
        for kl=1:size(x,1)
          ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
          set(ht,'HorizontalAlignment','right')
        end
        return
    end
    end

end

end