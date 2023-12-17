function lambdas=bisection(thetaz,l,a,b)

syms lambdaz

%     thetaz=lambdaz^2+2*lambdaz; %funcao a minimizar
%     l=0.01;
%     intervalo(1,1)=-3;intervalo(1,2)=6; %intervalo inicial de incerteza
%     a(1)=intervalo(1,1);b(1)=intervalo(1,2);

%% Passos principais

theta(1,1)=subs(thetaz,lambdaz,a(1));
theta(1,2)=subs(thetaz,lambdaz,b(1));

difftheta=diff(thetaz);

n=ceil(log(l/(b(1)-a(1)))/log(.5));

for k=1:n

    lambda(k)=.5*(a(k)+b(k));
    diff_lambda(k)=subs(difftheta,lambdaz,lambda(k));
   
    if diff_lambda(k)==0
        lambdas=(a(k+1)+b(k+1))/2;
        return
        
    elseif diff_lambda(k)>0 % Passo 2
        
        a(k+1)=a(k);
        b(k+1)=lambda(k);
        
        theta(k+1,1)=subs(thetaz,lambdaz,a(k+1));
        theta(k+1,2)=subs(thetaz,lambdaz,b(k+1));

    elseif diff_lambda(k)<0 % Passo 3
        
        a(k+1)=lambda(k);
        b(k+1)=b(k);
 
        theta(k+1,1)=subs(thetaz,lambdaz,a(k+1));
        theta(k+1,2)=subs(thetaz,lambdaz,b(k+1));

    end
    
    if isinf(diff_lambda(k)) || isnan(diff_lambda(k))
        
        error('Vicio no Bisection!')
        
    end

end

lambdas=(a(k+1)+b(k+1))/2;
end