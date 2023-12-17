function prob4c

format compact  
warning('off','all') 

syms x y

x1 = linspace(-.5,4,101);
x2 = linspace(-.5,4,101);

[X1, X2] = meshgrid(x1,x2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Funcao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = X1.^2-X1.*X2+2*X2.^2-4*X1-5*X2;
f = x^2-x*y+2*y^2-4*x-5*y;

G1 = X1+2*X2-6;
g1 = x+2*y-6;

G2 = X2-2;
g2 = y-2;

%% Solucao analitica
% Condicoes necessarias de primeira ordem
syms beta1 beta2

Af = f + beta1*g1 + beta2*g2;
dfdx1 = diff(Af,x);dfdx2 = diff(Af,y);  

x1b = 0;   x2b = 0;  beta1b = 0;   beta2b = 0;
fb = 0;   g1b = 0;   g2b = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CASO A : beta1 = 0, beta2 = 0
%%%%  Verificar g1 < 0 and g2 < 0 na solucao - CONSTRANGIMENTOS INACTIVOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n')
solCasea = solve(subs(dfdx1,{beta1,beta2},{0,0}),subs(dfdx2,{beta1,beta2},{0,0}),'x,y');

x1s = solCasea.x;
x2s = solCasea.y;
beta1s = subs(beta1,0);
beta2s = subs(beta2,0);

fs = subs(f,{x,y,beta1,beta2},{x1s,x2s,0,0});
g1s = subs(g1,{x,y,beta1,beta2},{x1s,x2s,0,0});
g2s = subs(g2,{x,y,beta1,beta2},{x1s,x2s,0,0});

fprintf('*******************************\n')
fprintf('As Solucoes para *** CASO A *** \n')
fprintf('*******************************\n')
fprintf('   x1*   x2*    u1*    u2*    f*   g1   g2\n')
SOL=double([x1s x2s beta1s beta2s fs g1s g2s]);
for i=1:length(SOL)-6;disp(SOL(1,i:i+6));end

%% Melhores solucoes para o CASO A
for i = 1:length(x1s)
    x1sv = double((x1s(i)));
    x2sv = double((x2s(i)));
    beta1v = double(beta1s(i));
    beta2v = double(beta2s(i));
    g1sv = double(g1s(i));
    g2sv = double(g2s(i));
    fsv = double(fs(i));
    if isreal(x1sv) && isreal(x2sv) && isreal(beta1v) && isreal(beta2v)
        if x1sv >= 0 && x2sv >= 0 
            if g1sv <= 0 && g2sv <= 0;if beta1v >= 0 && beta2v >= 0;if fsv < fb;x1b = x1sv;x2b = x2sv;beta1b = beta1v;beta2b = beta2v;fb = fsv;g1b = g1sv;g2b = g2sv;end;end;end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CASO B) : g1 = 0, beta2 = 0
%%%%  Verificar beta1 > 0 and g2 = 0 na solucao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n')
clear x1s x2s beta1s beta2s fs g1s g2s
solCaseb = solve(subs(dfdx1,{beta2},{0}),subs(dfdx2,{beta2},{0}),g1,'x,y,beta1');

x1s = solCaseb.x;
x2s = solCaseb.y;
beta1s = solCaseb.beta1;
beta2s = subs(beta2,0)*ones(size(beta1s),1);

fs = subs(f,{x,y,beta1,beta2},{x1s,x2s,beta1s,beta2s});
g1s = subs(g1,{x,y,beta1,beta2},{x1s,x2s,beta1s,beta2s});
g2s = subs(g2,{x,y,beta1,beta2},{x1s,x2s,beta1s,beta2s});

fprintf('*******************************\n')
fprintf('As Solucoes para *** CASO B *** \n')
fprintf('*******************************\n')
fprintf('     x1*     x2*        u1*       u2*      f*      g1        g2\n')
SOL=double([x1s x2s beta1s beta2s fs g1s g2s]);
for i=1:length(SOL)-6;disp(SOL(1,i:i+6));end

%% Melhores solucoes para o CASO B
for i = 1:length(x1s)
    x1sv = double((x1s(i)));
    x2sv = double((x2s(i)));
    beta1v = double(beta1s(i));
    beta2v = double(beta2s(i));
    g1sv = double(g1s(i));
    g2sv = double(g2s(i));
    fsv = double(fs(i));
    if isreal(x1sv) && isreal(x2sv) && isreal(beta1v) && isreal(beta2v)
        if x1sv >= 0 && x2sv >= 0 
            if g1sv <= 0 && g2sv <= 0;if beta1v >= 0 && beta2v >= 0;if fsv < fb;x1b = x1sv;x2b = x2sv;beta1b = beta1v;beta2b = beta2v;fb = fsv;g1b = g1sv;g2b = g2sv;end;end;end
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CASO C : beta1 = 0, g2 = 0
%%%%  Verificar g1 < 0 and beta2 > 0 na solucao 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n')
clear x1s x2s beta1s beta2s fs g1s g2s

solCasec = solve(subs(dfdx1,{beta1},{0}),subs(dfdx2,{beta1},{0}),g2,'x,y,beta2');

x1s = solCasec.x;
x2s = solCasec.y;
beta2s = solCasec.beta2;
beta1s = subs(beta1,0)*ones(size(beta2s),1);

fs = subs(f,{x,y,beta1,beta2},{x1s,x2s,beta1s,beta2s});
g1s = subs(g1,{x,y,beta1,beta2},{x1s,x2s,beta1s,beta2s});
g2s = subs(g2,{x,y,beta1,beta2},{x1s,x2s,beta1s,beta2s});

fprintf('*******************************\n')
fprintf('As Solucoes para *** CASO C *** \n')
fprintf('*******************************\n')
fprintf('   x1*     x2*  u1*    u2*  f*   g1    g2\n')
SOL=double([x1s x2s beta1s beta2s fs g1s g2s]);
for i=1:length(SOL)-6;disp(SOL(1,i:i+6));end

%% Melhores solucoes para o CASO C
for i = 1:length(x1s)
    x1sv = double((x1s(i)));
    x2sv = double((x2s(i)));
    beta1v = double(beta1s(i));
    beta2v = double(beta2s(i));
    g1sv = double(g1s(i));
    g2sv = double(g2s(i));
    fsv = double(fs(i));
    if isreal(x1sv) && isreal(x2sv) && isreal(beta1v) && isreal(beta2v)
        if x1sv >= 0 && x2sv >= 0 
            if g1sv <= 0 && g2sv <= 0;if beta1v >= 0 && beta2v >= 0;if fsv < fb;x1b = x1sv;x2b = x2sv;beta1b = beta1v;beta2b = beta2v;fb = fsv;g1b = g1sv;g2b = g2sv;end;end;end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CASO D : g1 = 0, g2 = 0
%%%%  Verificar beta1 > 0 and beta2 > 0 na solucao 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n')
clear x1s x2s beta1s beta2s fs g1s g2s

solCased = solve(subs(dfdx1),subs(dfdx2),g1,g2,'x,y,beta1,beta2');

x1s = solCased.x;
x2s = solCased.y;
beta1s = solCased.beta1;
beta2s = solCased.beta2;

fs = subs(f,{x,y,beta1,beta2},{x1s,x2s,beta1s,beta2s});
g1s = subs(g1,{x,y,beta1,beta2},{x1s,x2s,beta1s,beta2s});
g2s = subs(g2,{x,y,beta1,beta2},{x1s,x2s,beta1s,beta2s});

fprintf('*******************************\n')
fprintf('As Solucoes para *** CASO D *** \n')
fprintf('*******************************\n')
fprintf('  x1*    x2*   u1*    u2*    f*   g1   g2\n')
SOL=double([x1s x2s beta1s beta2s fs g1s g2s]);
for i=1:length(SOL)-6;disp(SOL(1,i:i+6));end

%% Melhores solucoes para o CASO D
for i = 1:length(x1s)
    x1sv = double((x1s(i)));
    x2sv = double((x2s(i)));
    beta1v = double(beta1s(i));
    beta2v = double(beta2s(i));
    g1sv = double(g1s(i));
    g2sv = double(g2s(i));
    fsv = double(fs(i));
    if isreal(x1sv) && isreal(x2sv) && isreal(beta1v) && isreal(beta2v)
        if x1sv >= 0 && x2sv >= 0 
            if g1sv <= 0 && g2sv <= 0;if beta1v >= 0 && beta2v >= 0;if fsv < fb;x1b = x1sv;x2b = x2sv;beta1b = beta1v;beta2b = beta2v;fb = fsv;g1b = g1sv;g2b = g2sv;end;end;end
        end
    end
end

%%%%  Coleccionar as melhores solucoes para os CASOS A, B, C e D
for i = 1:length(x1s)
    x1sv = double((x1s(i)));
    x2sv = double((x2s(i)));
    beta1v = double(beta1s(i));
    beta2v = double(beta2s(i));
    g1sv = double(g1s(i));
    g2sv = double(g2s(i));
    fsv = double(fs(i));
    if isreal(x1sv) && isreal(x2sv) && isreal(beta1v) && isreal(beta2v)
        if x1sv >= 0 && x2sv >= 0 
            if g1sv <= 0 && g2sv <= 0;if beta1v >= 0 && beta2v >= 0;if fsv < fb;x1b = x1sv;x2b = x2sv;beta1b = beta1v;beta2b = beta2v;fb = fsv;g1b = g1sv;g2b = g2sv;end;end;end
        end
    end
end

fprintf('*******************************\n')
fprintf('         MELHOR SOLUCAO:       \n')
fprintf('*******************************\n')
fprintf('      x1*      x2*        beta1*      beta2*    f*       g1        g2\n'), ...
   disp(double([x1b x2b beta1b beta2b fb g1b g2b]))

syms lambdaz 
syms x1 x2 real

thetas=x1^2-x1*x2+2*x2^2-4*x1-5*x2; %funcao a minimizar
thetas_desig=[x1+2*x2-6;x2-2];ndes=size(thetas_desig,1); %sujeito a desigualdade
thetas_ig=0;ni=size(thetas_ig,1); %sujeito a igualdade

gradfx=sym(zeros(2,1));

% Gradiente
dfdx1=diff(thetas,x1);dfdx2=diff(thetas,x2);
gradfx(1,1)=dfdx1;gradfx(2,1)=dfdx2;

%gradientes das igualdades/desigualdades

gradg=sym(zeros(2,ndes));
for j=1:ndes
gradg(:,j)=[diff(thetas_desig(j),x1);diff(thetas_desig(j),x2)];
end

gradi=sym(zeros(2,ni));
for j=1:ni
gradi(:,j)=[diff(thetas_ig(j),x1);diff(thetas_ig(j),x2)];
end

%escolher os multiplicadores de lagrange do problema
syms u1 v1 u2 v2 u3 v3

u=[u1,u2]; nu=size(u,2); v=[v1]; nv=size(v,2);

sumg=0;sumg1=sym(zeros(nu,1));
    for i=1:nu
    sumg=sumg+u(i)*gradg(:,i)';
    sumg1(i,1)=u(i)*thetas_desig(i,1);
    end
    somaG=vpa(sumg);
    sumg1=vpa(sumg1);
    
sumh=0;
    for i=1:nv
    sumg=sumh+v(i)*gradi(:,i)';
    end
    somaH=vpa(sumg);
    
eq1=gradfx(1,1)+somaG(1)+somaH(1);
eq2=gradfx(2,1)+somaG(2)+somaH(2);

fprintf('\n')
fprintf('-----------------------------------------------------------------------------\n')
fprintf('Condicoes necessarias de KKT:\n%s=0\n%s=0\n',char(eq1),char(eq2));
fprintf('%s=%g\n%s=%g\n',char(sumg1(1,1)),0,char(sumg1(2,1)),0);
for ii=1:nu
fprintf('%s>=%g\n',char(u(ii)),0);
end
fprintf('-----------------------------------------------------------------------------\n')
fprintf('\n')

ver_ponto=[x1b x2b];%ver_ponto=[2.6 1.7];

[su1,su2]=solve(subs(eq1,[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]),subs(eq2,[x1,x2],[ver_ponto(1,1),ver_ponto(1,2)]),u1,u2);
su1=double(su1);su2=double(su2);

kktcnd31=subs(sumg1(1,1),[x1,x2,u1,u2],[ver_ponto(1,1),ver_ponto(1,2),su1,su2]);
kktcnd32=subs(sumg1(2,1),[x1,x2,u1,u2],[ver_ponto(1,1),ver_ponto(1,2),su1,su2]);
if ~strcmp(class(kktcnd31),'sym');kktcnd31=double(kktcnd31);end
if ~strcmp(class(kktcnd32),'sym');kktcnd32=double(kktcnd32);end

if (su1>=0 & ~isempty(su1)) & ~isempty(su2) & kktcnd31==0 & kktcnd32==0
    fprintf('-----------------------------------------------------------------------------\n')
    fprintf('O ponto (%6.3f,%6.3f) satisfaz as condicoes de KKT!\n',x1b,x2b);
    fprintf('u1=%6.3f;u2=%6.3f\n',su1,su2);
    fprintf('-----------------------------------------------------------------------------\n')
else
    fprintf('-----------------------------------------------------------------------------\n')
    fprintf('O ponto (%6.3f,%6.3f) nao satisfaz as condicoes de KKT!\n',x1b,x2b);
    fprintf('-----------------------------------------------------------------------------\n')
end

end