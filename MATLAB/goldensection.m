function lambdas=goldensection(thetaz,l,a,b,alpha,lambda,miu)

syms lambdaz

%% Passos principais

theta(1,1)=subs(thetaz,lambdaz,lambda(1));
theta(1,2)=subs(thetaz,lambdaz,miu(1));

% figure(K)
% ezplot(thetaz,[a(1) b(1)]);hold on; % Plot da funcao no intervalo inicial

k=1;
while b(k)-a(k)>l
    
    if theta(k,1)>theta(k,2) % Passo 2
        
        a(k+1)=lambda(k);
        b(k+1)=b(k);
        lambda(k+1)=miu(k);
        miu(k+1)=a(k+1)+alpha*(b(k+1)-a(k+1));
        
        theta(k+1,1)=theta(k,2);
        theta(k+1,2)=subs(thetaz,lambdaz,miu(k+1));
         
%         fprintf('-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n')
%         fprintf('a(%6.0f)= %6.4f \n',k+1,a(k+1))
%         fprintf('b(%6.0f)= %6.4f \n',k+1,b(k+1))
%         fprintf('miu(%6.0f)= %6.4f \n',k+1,miu(k+1))
%         fprintf('theta1(%6.0f)= %6.4f \n',k+1,theta(k+1,1))
%         fprintf('theta2(%6.0f)= %6.4f \n',k+1,theta(k+1,2))

        
    else % Passo 3
        
        a(k+1)=a(k);
        b(k+1)=miu(k);
        miu(k+1)=lambda(k);
        lambda(k+1)=a(k+1)+(1-alpha)*(b(k+1)-a(k+1));
        
        theta(k+1,1)=subs(thetaz,lambdaz,lambda(k+1));
        theta(k+1,2)=theta(k,1);
        
%         fprintf('-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n')
%         fprintf('a(%6.0f)= %6.4f \n',k+1,a(k+1))
%         fprintf('b(%6.0f)= %6.4f \n',k+1,b(k+1))
%         fprintf('lambdaGS(%6.0f)= %6.4f \n',k+1,lambda(k+1))
%         fprintf('theta1(%6.0f)= %6.4f \n',k+1,theta(k+1,1))
%         fprintf('theta2(%6.0f)= %6.4f \n',k+1,theta(k+1,2))
 
    end

    k=k+1; % Passo 4

end

% fprintf('-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n')

lambdas=(lambda(k)+miu(k))/2;

end