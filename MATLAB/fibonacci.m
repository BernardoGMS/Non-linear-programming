function lambdas=fibonacci(thetaz,l,a,b,alpha,lambda,miu,epsilon)
clc
syms lambdaz

    epsilon=0.01;
    thetaz=lambdaz^2+2*lambdaz; %funcao a minimizar
    l=0.2;
    intervalo(1,1)=-3;intervalo(1,2)=5; %intervalo inicial de incerteza
    a(1)=intervalo(1,1);b(1)=intervalo(1,2);

%% Passos principais

FIB(1)=0;FIB(2)=1;for i=1:100;FIB(i+2)=FIB(i)+FIB(i+1);end;
nn=(b(1)-a(1))/l;j=1;

while FIB(j)<nn
    n=j-1;j=j+1;
end

lambda(1)=a(1)+(FIB(n-1)/FIB(n+1))*(b(1)-a(1));
miu(1)=a(1)+(FIB(n)/FIB(n+1))*(b(1)-a(1));

theta(1,1)=subs(thetaz,lambdaz,lambda(1));
theta(1,2)=subs(thetaz,lambdaz,miu(1));

k=1;
while k<n

    disp(a)
    disp(b)
    disp(lambda)
    disp(miu)
    disp(theta(:,1))
    disp(theta(:,2))
    if theta(k,1)>theta(k,2) % Passo 2
        
        a(k+1)=lambda(k);
        b(k+1)=b(k);
        lambda(k+1)=miu(k);
        miu(k+1)=a(k+1)+(FIB(n-k)/FIB(n-k+1))*(b(k+1)-a(k+1));
        
        if k==n-2; % Passo 5
            
            lambda(n)=lambda(n-1);
            miu(n)=lambda(n-1)+epsilon;
            
            if subs(thetaz,lambdaz,lambda(n))>subs(thetaz,lambdaz,miu(n))
                
                a(n)=lambda(n);
                b(n)=b(n-1);
                lambdas=(a(n)+b(n))/2;
                
                return
                
            elseif subs(thetaz,lambdaz,lambda(n))<=subs(thetaz,lambdaz,miu(n))
                
                a(n)=a(n-1);
                b(n)=lambda(n);
                lambdas=(a(n)+b(n))/2;
                
                return
                
            end
            
            
        else
            
            theta(k+1,1)=subs(thetaz,lambdaz,lambda(k+1));
            theta(k+1,2)=subs(thetaz,lambdaz,miu(k+1));
            k=k+1;
            
        end
        
    else % Passo 3
        
        a(k+1)=a(k);
        b(k+1)=miu(k);
        miu(k+1)=lambda(k);
        lambda(k+1)=a(k+1)+(FIB(n-k-1)/FIB(n-k+1))*(b(k+1)-a(k+1));
        
        if k==n-2; % Passo 5
            
            lambda(n)=lambda(n-1);
            miu(n)=lambda(n-1)+epsilon;
            
            if subs(thetaz,lambdaz,lambda(n))>subs(thetaz,lambdaz,miu(n))
                
                a(n)=lambda(n);
                b(n)=b(n-1);
                lambdas=(a(n)+b(n))/2
                k
                return
                
            elseif subs(thetaz,lambdaz,lambda(n))<=subs(thetaz,lambdaz,miu(n))
                
                a(n)=a(n-1);
                b(n)=lambda(n);
                lambdas=(a(n)+b(n))/2
                k
                return
                
            end
            
            
        else
            
            theta(k+1,1)=subs(thetaz,lambdaz,lambda(k+1));
            theta(k+1,2)=subs(thetaz,lambdaz,miu(k+1));
            k=k+1;
            
        end
 
    end

end


end