function [x,VIOLG,VIOLI,miu,miv,v,u,theta,THETA]=prob3c(metodo)
syms x1 x2

% Inicializacao
[x,miu,miv,v,u,thetaz,igthetaz,desigthetaz,alpha]=prob3b;
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
elseif strcmp(metodo,'DFP')
[xx,THETA(k,1)]=DFP(alpha);
else
error('Metodo desconhecido!!!!')
end

for i=1:25

x(k+1,:)=xx(end,:);
theta(k+1,1)=subs(thetaz,[x1,x2],[x(k+1,1),x(k+1,2)]);

plot(x(k:k+1,1),x(k:k+1,2),'ro-');hold on;

VIOLG(k+1,:)=max(subs(desigthetaz,[x1,x2],[x(k+1,1),x(k+1,2)]),0);
VIOLI(k+1,:)=max(abs(subs(igthetaz,[x1,x2],[x(k+1,1),x(k+1,2)])));
activeG(1,:)=zeros(1,Ndes);activeI(1,:)=zeros(1,Ndes);

%% Verificacao das violacoes das desigualdades

for j=1:Ndes
if VIOLG(k+1,j)>epsilon && activeG(k,j)==0
    activeG(k+1,j)=1;
    u(k+1,j)=u(k,j)+max((2*miu(k,j))*subs(desigthetaz(j,1),[x1,x2],[x(k+1,1),x(k+1,2)]),0);
    miu(k+1,j)=10*miu(1,j);fprintf('Activa %6.2f miu(%s,%s)...    ',miu(k+1,j),int2str(k+1),int2str(j))
elseif VIOLG(k+1,j)>epsilon && activeG(k,j)==1
    activeG(k+1,j)=1;
    u(k+1,j)=u(k,j)+max((2*miu(k,j))*subs(desigthetaz(j,1),[x1,x2],[x(k+1,1),x(k+1,2)]),0);
    miu(k+1,j)=10*miu(k,j);fprintf('Actualiza %6.2f miu(%s,%s)...    ',miu(k+1,j),int2str(k+1),int2str(j)) 
else activeG(k+1,j)=0;miu(k+1,j)=miu(k,j);u(k+1,j)=u(k,j);fprintf('Nao actualiza %6.2f miu(%s,%s)...    ',miu(k+1,j),int2str(k+1),int2str(j))
end
fprintf('\n')
end

%% Verificacao das violacoes das igualdades

for j=1:Nig
if VIOLI(k+1,j)>epsilon && activeI(k,j)==0
    activeI(k+1,j)=1;v(k+1,j)=v(k,j)+2*miv(k,j)*subs(igthetaz(j,1),[x1,x2],[x(k,1),x(k,2)]);
    miv(k+1,j)=10*miv(1,j);fprintf('Activa %6.2f miu(%s,%s)...    ',miv(k+1,j),int2str(k+1),int2str(j))
elseif VIOLI(k+1,j)>epsilon && activeI(k,j)==1
    activeI(k+1,j)=1;v(k+1,j)=v(k,j)+2*miv(k,j)*subs(igthetaz(j,1),[x1,x2],[x(k,1),x(k,2)]);
    miv(k+1,j)=10*miv(k,j);fprintf('Actualiza %6.2f miu(%s,%s)...    ',miv(k+1,j),int2str(k+1),int2str(j)) 
else activeI(k+1,j)=0;miv(k+1,j)=miv(k,j);v(k+1,j)=v(k,j);fprintf('Nao actualiza %6.2f miu(%s,%s)...    ',miv(k+1,j),int2str(k+1),int2str(j))
end
fprintf('\n')
end

    k=k+1;
    
    % Calcular os termos das desigualdades
    somaG=0;
    for j=1:Ndes
      somaG=somaG+miu(k,j)*(desigthetaz(j,1)+max(u(k,j)/(2*miu(k,j)),0))^2;
    end
    
    alpha=thetaz+somaG-sum(u(k,:).^2./(4*miu(k,:)))+v(k,:)*igthetaz+miv(k,:)*igthetaz.^2;
    fprintf('Alpha de k=%s...  \n',int2str(k))
    
    kk=0;
    for j=1:Ndes
        fprintf('VIOLG(k,j) %E k=%s...  \n',VIOLG(k,j),int2str(k))
        if VIOLG(k,j)<epsilon
        kk=kk+1;
        end
    end
    
    for j=1:Nig
        fprintf('VIOLI(k,j) %E k=%s...  \n',VIOLI(k,j),int2str(k))
        if VIOLI(k,j)<epsilon
        kk=kk+1;
        end
    end
    
    fprintf('\n')
    fprintf('FUNCAO DESCONSTRANGIDA F(xk,vi,ui):\n')
    pretty(vpa(alpha))
    
%% Procurar nova solucao
    if strcmp(metodo,'GC')
    [xx,THETA(k,1)]=gradientes_conj_FleRee(alpha);
    elseif strcmp(metodo,'DFP')
    [xx,THETA(k,1)]=DFP(alpha);
    else
    error('Metodo desconhecido!!!!')
    end
    
    if norm(x(k,:)-x(k-1,:))<epsilon && k>2 && kk==Nig+Ndes
        for kl=1:size(x,1)
          ht = text(x(kl,1),x(kl,2),['x',int2str(kl)]);
          set(ht,'HorizontalAlignment','right','FontSize',14)
        end
        return
    end
    
end

end