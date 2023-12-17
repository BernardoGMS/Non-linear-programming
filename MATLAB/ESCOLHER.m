clc
clear all
close all

fprintf('______________________________________________________________\n')
fprintf('       METODOS COMPUTACIONAIS E DE OTIMIZACAO 2014/2015\n')
fprintf('          Goncalo Moniz Silva Bernardo (ist90893)\n')
fprintf('______________________________________________________________\n')
fprintf('\n')

nr_problema=2; %Escolha do problema

problema={'problema1';'problema2';'problema3';'problema4';'problema5'};
escolhe_problema=char(problema(nr_problema));

switch lower(escolhe_problema)
    case {'problema1'}
        
        fprintf('______________________________________________________________\n')
        fprintf('                        EXERCICIO 1\n')
        fprintf('______________________________________________________________\n') 

        [intervalo,lambda,miu,theta,n,compincert]=prob1;
        
        fprintf('--------------------------------------------------------------- \n')
        fprintf('Solucao optima encontrada no intervalo [%6.2f %6.2f]: %6.5f \n',intervalo(end,1),intervalo(end,2),(theta(end,1)+theta(end,2))/2)
        fprintf('Incerteza final e numero minimo de observacoes: %6.3f/%6.0f \n',compincert(end),n)
        
        clearvars -except intervalo lambda miu theta compincert
        
    case 'problema2'
        
        alinea='d';
        
        if strcmp(alinea,'a')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 2, ALINEA a\n')
        fprintf('                 Condicoes Necessarias de Minimo\n')
        fprintf('______________________________________________________________\n')
        prob2a
        elseif strcmp(alinea,'b')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 2, ALINEA b\n')
        fprintf('                    METODO DE HOOKE-JEEVES\n')
        fprintf('______________________________________________________________\n')
        [x,y,D,lambda_,theta,delta,ymais,ymenos]=prob2b;
        fprintf('--------------------------------------------------------------- \n')
        fprintf('SOLUCAO OPTIMA ENCONTRADA, PELO METODO DE HOOKE-JEEVES:\n       x=[%6.3f %6.3f]; f=%6.4f \n',x(end,1),x(end,2),theta(end))
        fprintf('--------------------------------------------------------------- \n')
        elseif strcmp(alinea,'c')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 2, ALINEA c\n')
        fprintf('                    METODO DE ROSENBROCK\n')
        fprintf('______________________________________________________________\n')
        [x,y,d,lambda_,THETA,delta,ymais,ymenos]=prob2c;
        fprintf('--------------------------------------------------------------- \n')
        fprintf('SOLUCAO OPTIMA ENCONTRADA, PELO METODO DE ROSENBROCK:\n       x=[%6.3f %6.3f]; f=%6.4f \n',x(end,1),x(end,2),THETA(end))
        fprintf('--------------------------------------------------------------- \n')
        elseif strcmp(alinea,'d')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 2, ALINEA d\n')
        fprintf('              METODO DOS GRADIENTES CONJUGADOS\n')
        fprintf('______________________________________________________________\n')    
        [x,y,d,THETA,alpha]=prob2d;
        fprintf('--------------------------------------------------------------- \n')
        fprintf('SOLUCAO OPTIMA ENCONTRADA, PELO METODO DOS GRADIENTES CONJUGADOS:\n       x=[%6.3f %6.3f]; f=%6.4f \n',x(end,1),x(end,2),THETA(end))
        fprintf('--------------------------------------------------------------- \n')
        elseif strcmp(alinea,'e')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 2, ALINEA e\n')
        fprintf('                        METODO DE DFP\n')
        fprintf('______________________________________________________________\n') 
        [x,y,d,THETA,D,lambdaj]=prob2e;
        fprintf('--------------------------------------------------------------- \n')
        fprintf('SOLUCAO OPTIMA ENCONTRADA, PELO METODO DE DFP:\n       x=[%6.3f %6.3f]; f=%6.4f \n',x(end,1),x(end,2),THETA(end))
        fprintf('--------------------------------------------------------------- \n')
        else
        fprintf('\n')
        fprintf('Alinea inexistente!\n')
        fprintf('\n')
        end
        
        clearvars -except x y d lambda_ theta delta ymais ymenos THETA alpha D lambdaj
        
    case 'problema3'
        
        alinea='c';
        
        if strcmp(alinea,'a')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 3, ALINEA a\n')
        fprintf('                  Condicoes Necessarias de KKT\n')
        fprintf('______________________________________________________________\n')
        prob3a
        elseif strcmp(alinea,'b')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 3, ALINEA b\n')
        fprintf('               Formulacao Langrangeana Aumentada\n')
        fprintf('______________________________________________________________\n')
        [x,miu,miv,v,u,thetaz,igthetaz,desigthetaz,alpha]=prob3b;
        elseif strcmp(alinea,'c')
        metodo='GC';
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 3, ALINEA c\n')
        fprintf('        Metodo dos Gradientes Conjugados de Fletcher-Reeves\n')
        fprintf('______________________________________________________________\n')
        fprintf('\n')
        [x,VIOLG,VIOLI,miu,miv,v,u,theta,THETA]=prob3c(metodo);
        fprintf('--------------------------------------------------------------- \n')
        fprintf('SOLUCAO OPTIMA ENCONTRADA, PELO METODO DOS GC:\n       x=[%6.3f %6.3f]; f=%6.4f \n',x(end,1),x(end,2),theta(end))
        fprintf('--------------------------------------------------------------- \n')
%         prob3e(x(end,:));
        elseif strcmp(alinea,'d')
        metodo='DFP';
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 3, ALINEA d\n')
        fprintf('                        Metodo de DFP\n')
        fprintf('______________________________________________________________\n')
        fprintf('\n')
        [x,VIOLG,VIOLI,miu,miv,v,u,theta,THETA]=prob3c(metodo);
        fprintf('--------------------------------------------------------------- \n')
        fprintf('SOLUCAO OPTIMA ENCONTRADA, PELO METODO DE DPF:\n       x=[%6.3f %6.3f]; f=%6.4f \n',x(end,1),x(end,2),THETA(end))
        fprintf('--------------------------------------------------------------- \n')
%         prob3e(x(end,:));
        else
        fprintf('\n')
        fprintf('Alinea inexistente!')
        fprintf('\n')
        end
         
        clearvars -except x mi v f alphas theta VIOLG VIOLI miu miv u
        
    case 'problema4'
        
        alinea='a';
        
        if strcmp(alinea,'a')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 4, ALINEA a\n')
        fprintf('                  Resolucao Grafica do Problema\n')
        fprintf('______________________________________________________________\n')
        prob4a
        elseif strcmp(alinea,'b')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 4, ALINEA b\n')
        fprintf('                  Condicoes Necessarias de KKT\n')
        fprintf('______________________________________________________________\n')
        prob4b;
        elseif strcmp(alinea,'c')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 4, ALINEA c\n')
        fprintf('           Verificacao das Condicoes Necessarias de KKT\n')
        fprintf('______________________________________________________________\n')
        prob4c;
        elseif strcmp(alinea,'d') 
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 4, ALINEA d\n')
        fprintf('              Formulacao de Penalidade do Problema\n')
        fprintf('______________________________________________________________\n')
        [x,miu,miv,thetaz,igthetaz,desigthetaz,alpha]=prob4d;
        elseif strcmp(alinea,'e')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 4, ALINEA d\n')
        fprintf('                     Metodo de ROSENBROCK\n')
        fprintf('______________________________________________________________\n')
        metodo='RS';
        [x,VIOLG,VIOLI,miu,miv,theta,THETA]=prob4f(metodo);
        fprintf('--------------------------------------------------------------- \n')
        fprintf('SOLUCAO OPTIMA ENCONTRADA, PELO METODO DE RS:\n       x=[%6.3f %6.3f]; f=%6.4f \n',x(end,1),x(end,2),theta(end))
        fprintf('--------------------------------------------------------------- \n')
        elseif strcmp(alinea,'f')
        fprintf('______________________________________________________________\n')
        fprintf('                    EXERCICIO 4, ALINEA f\n')
        fprintf('                Metodo dos GRADIENTES CONJUGADOS\n')
        fprintf('______________________________________________________________\n')
        metodo='GC';
        [x,VIOLG,VIOLI,miu,miv,theta,THETA]=prob4f(metodo);
        fprintf('--------------------------------------------------------------- \n')
        fprintf('SOLUCAO OPTIMA ENCONTRADA, PELO METODO DOS GC:\n       x=[%6.3f %6.3f]; f=%6.4f \n',x(end,1),x(end,2),theta(end))
        fprintf('--------------------------------------------------------------- \n')
        else
        fprintf('\n')
        fprintf('Alinea inexistente!')
        fprintf('\n')
        end
 
        clearvars -except x f alphas theta VIOLG VIOLI miu miv THETA
        
    case 'problema5'
        
        alinea='c';
        
        if strcmp(alinea,'ab') 
        fprintf('______________________________________________________________\n')
        fprintf('                   EXERCICIO 5, ALINEAS a E b\n')
        fprintf('              CONDICOES DE 1� E 2� ORDEM DE MINIMO\n')
        fprintf('______________________________________________________________\n')
        thetas=prob5ab;  
        elseif strcmp(alinea,'c')
        fprintf('______________________________________________________________\n')
        fprintf('                   EXERCICIO 5, ALINEA c\n')
        fprintf('                        Metodo de DFP\n')
        fprintf('______________________________________________________________\n')
        [x,y,d,THETA,D,lambdaj]=prob5c;
        end
        
        clearvars -except x y d D lambdaj THETA
        
    otherwise
        
    error('Problema desconhecido!');
        
end

fprintf('=======================================================================\n')
fprintf('                       %s\n',datestr(now))
fprintf('=======================================================================\n')