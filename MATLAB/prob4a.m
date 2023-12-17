function prob4a

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

x1u = 2;  x2u = 2;
fval = subs(f,{x,y},{x1u, x2u});
cvalues = linspace(-11,fval,0);
[C1,h1] = contour(x1,x2,F,[cvalues -10.78],'g--');clabel(C1,h1);

grid on
title('Problema 4a: RESOLUCAO GRAFICA')
xlabel('x_1')
ylabel('x_2')

hold on
[C1,h1] = contour(x1,x2,G1,[0 0],'r-');
clabel(C1,h1);set(h1,'LineWidth',2)

[C2,h2] = contour(x1,x2,G2,[0 0],'b-');
clabel(C2,h2);set(h2,'LineWidth',2);

plot(x1b,x2b,'ro','MarkerSize',15, ...
    'MarkerFaceColor','y')

ht = text(2.6,1.7,['x=(',num2str(2.6),',',num2str(1.7),')']);set(ht,'HorizontalAlignment','right')

h3 = line([0 4],[0 0],'Color','k', ...
    'LineWidth',2,'LineStyle','-');
h4 = line([0 0],[0 4],'Color','k', ...
    'LineWidth',2,'LineStyle','-');
axis square

%% Gradiente da Funcao em xoptimo
dfdx1 = diff(f,x);
dfdx2 = diff(f,y);
grad =[dfdx1;dfdx2];
% gradval = double(subs(grad,{x,y},{x1b,x2b}));
gradval = double(subs(grad,{x,y},{2.6,1.7}));
delx1g = [0 -0.2];
delx2_g =  (gradval(2,1)/gradval(1,1))*delx1g;
% x1_g = x1b + delx1g;
% x2_g = x2b + delx2_g;
x1_g = 2.6 + delx1g;
x2_g = 1.7 + delx2_g;

hg1 = line(x1_g, x2_g,'Color','g', ...
     'LineWidth',1.5,'LineStyle','>');
hg2 = line(x1_g, x2_g,'Color','g', ...
     'LineWidth',1.5,'LineStyle','--');

 %% Gradiente de G1 em xoptimo
dg1dx1 = diff(g1,x);
dg1dx2 = diff(g1,y);
gradh =[dg1dx1;dg1dx2];
gradhval = double(subs(gradh,{x,y},{x1b,x2b}));
delx1g_h = [0 +0.15];
delx2_gh =  (gradhval(2,1)/gradhval(1,1))*delx1g_h;
x1_gh = x1b + delx1g_h;
x2_gh = x2b + delx2_gh;

hg3 = line(x1_gh, x2_gh,'Color','m', ...
     'LineWidth',2,'LineStyle','<');
hg4 = line(x1_gh, x2_gh,'Color','m', ...
     'LineWidth',2,'LineStyle','--');
 
%% Gradiente de G2 em xoptimo
dg2dx1 = diff(g2,x);
dg2dx2 = diff(g2,y);
gradh =[dg2dx1;dg2dx2];
gradhval = double(subs(gradh,{x,y},{x1b,x2b}));
delx1g_h = [0 +0.25];
delx2_gh =  (gradhval(2,1)/gradhval(1,1))*delx1g_h;
x1_gh = x1b + delx1g_h;
x2_gh = x2b + delx2_gh;

hg5 = line(x1_gh, x2_gh,'Color','r', ...
     'LineWidth',2,'LineStyle','>');
hg6 = line(x1_gh, x2_gh,'Color','r', ...
     'LineWidth',2,'LineStyle','--');
 
 legend('F(x)=x_1^2+x_1x_2+2x_2^2-4x_1-5x_2','x_1+2x_2-6<=0','x_2-2<=0',['(',num2str(2.6),',',num2str(1.7),')'],'x_1>=0','x_2>=0','Grad(F(x_o_p_t))','','Grad(G1(X_o_p_t))','','Grad(G2(X_o_p_t))','','Location','southeastoutside')

hold off

end
