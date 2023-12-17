function [x,VIOLG,VIOLI,miu,miv,v,u,theta,THETA]=prob4e(metodo)

syms x1 x2

% Inicializacao
[x,~,~,miu,miv,~,thetaz,igthetaz,desigthetaz,alpha]=prob4d;
Nig=size(igthetaz,1);Ndes=size(desigthetaz,1);

theta(1,1)=subs(thetaz,[x1,x2],[x(1,1),x(1,2)]);

epsilon=0.01;omega=.5;
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

for i=1:inf

x(k+1,:)=xx(end,:);
theta(k+1,1)=subs(thetaz,[x1,x2],[x(k+1,1),x(k+1,2)]);

VIOLG(k+1,:)=max(subs(desigthetaz,[x1,x2],[x(k+1,1),x(k+1,2)]),0);
VIOLI(k+1,:)=max(abs(subs(igthetaz,[x1,x2],[x(k+1,1),x(k+1,2)])));

%% Verificacao das violacoes das desigualdades
for j=1:Ndes
    
if VIOLG(k+1,j)>=VIOLG(k,j)
%     Rotina externa - Actualizacao do Multiplicador de Lagrange
%     u(k+1,j)=u(k,j)+max(2*miu(k,j)*subs(desigthetaz(j,1),[x1,x2],[x(k+1,1),x(k+1,2)]),-u(k,j));
u(k+1,j)=2*miu(k,j)*max(u(k,j)/(2*miu(k,j))+subs(desigthetaz(j,1),[x1,x2],[x(k+1,1),x(k+1,2)]),0);
% u(k+1,j)=max(u(k,j)+2*miu(k,j)*subs(desigthetaz(j,1),[x1,x2],[x(k,1),x(k,2)]),0);
    fprintf('Actualiza %6.2f u(%s,%s)...    ',u(k+1,j),int2str(k+1),int2str(j))
    fprintf('G %6.5f ...    ',subs(desigthetaz(j,1),[x1,x2],[x(k+1,1),x(k+1,2)]))
%     u(k+1,j)=u(k,j);
else u(k+1,j)=u(k,j); fprintf('Nao actualiza %6.2f u(%s,%s)...    ',u(k+1,j),int2str(k+1),int2str(j))
end

if norm(u(k+1,j)-u(k,j))<omega
    miu(k+1,j)=10*miu(k,j);fprintf('Actualiza %6.2f miu(%s,%s)...    ',miu(k+1,j),int2str(k+1),int2str(j))
else miu(k+1,j)=miu(k,j);fprintf('Nao actualiza %6.2f miu(%s,%s)...    ',miu(k+1,j),int2str(k+1),int2str(j))
end
fprintf('\n')
end

%% Verificacao das violacoes das igualdades

for j=1:Nig
    
if VIOLI(k+1,j)>=VIOLI(k,j)
    % Rotina externa - Actualizacao do Multiplicador de Lagrange
%     v(k+1,j)=v(k,j)+2*miv(k,j)*subs(igthetaz(j,1),[x1,x2],[x(k+1,1),x(k+1,2)]);
    v(k+1,j)=v(k,j)+2*miv(k,j)*subs(igthetaz(j,1),[x1,x2],[x(k,1),x(k,2)]);
    fprintf('Actualiza %6.2f v(%s,%s)...    ',v(k+1,j),int2str(k+1),int2str(j))
else v(k+1,j)=v(k,j);fprintf('Nao actualiza %6.2f v(%s,%s)...    ',v(k+1,j),int2str(k+1),int2str(j))
end

if norm(v(k+1,j)-v(k,j))<omega
    miv(k+1,j)=10*miv(k,j);fprintf('Actualiza %6.2f miv(%s,%s)...    ',miv(k+1,j),int2str(k+1),int2str(j))
else miv(k+1,j)=miv(k,j);fprintf('Nao actualiza %6.2f miv(%s,%s)...    ',miv(k+1,j),int2str(k+1),int2str(j))
end
fprintf('\n')
end

    k=k+1;
    
    % Calcular os termos das desigualdades
    somaG=0;
    for j=1:Ndes
%     somaG=somaG+miu(k,j)*max(subs(desigthetaz(j,1),[x1,x2],[x(k,1),x(k,2)])+u(k,j)/(2*miu(k,j)),0)^2;
      somaG=somaG+miu(k,j)*(desigthetaz(j,1)+max(u(k,j)/(2*miu(k,j)),0))^2;
    end
    
    alpha=thetaz+somaG-sum(u(k,:).^2./(4*miu(k,:)))+v(k,:)*igthetaz+miv(k,:)*igthetaz.^2;
    fprintf('Alpha de k=%s...  \n',int2str(k))
%     pretty(alpha)
    
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
    
    pretty(alpha)
    
%% Procurar nova solucao
    if strcmp(metodo,'GC')
    [xx,THETA(k,1)]=gradientes_conj_FleRee(alpha);
    elseif strcmp(metodo,'DFP')
    [xx,THETA(k,1)]=DFP(alpha);
    else
    error('Metodo desconhecido!!!!')
    end
    
    if norm(THETA(k,:)-THETA(k-1,:))<epsilon && k>2 && kk==Nig+Ndes
        return
    end

if i==10
    disp(i)
    return
end
end

end