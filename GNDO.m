
function [cgcurve,bestFitness,bestSol]=GNDO(obj,n,d,lb,ub,t)

%obj--------objective function
%c-------population size
%d-------dimension of problem
%lb-----the lower limit of the variables
%ub-----the upper limit of the variables
%t------the maximum number of function evaluations
%cgcurve---the record of the convergence curves
%bestobj--the optimal fitness value
%bestsol-------the optimal solution


% Initialise the population
for i=1:n
    x(i,:)=lb+(ub-lb).*rand(1,d); % Eq. 26
end

bestFitness = inf;

for it=1: 1 : t
    
    for i=1:n
        fitness(i) = obj(x(i,:));
        
        if fitness(i) < bestFitness
            bestSol = x(i,:);
        end
    end
    cgcurve(it)=bestFitness;
    
    
    mo= mean(x);
 

    for i=1:n
        a=randperm(n,1);
        b=randperm(n,1);
        c=randperm(n,1);
        while a==i | a==b | c==b | c==a |c==i |b==i
            a=randperm(n,1);
            b=randperm(n,1);
            c=randperm(n,1);
        end
        
        if fitness(a)<fitness(i)  %Eq. 24
            v1=x(a,:)-x(i,:);  
        else
            v1=x(i,:)-x(a,:);
        end
        
        if fitness(b)<fitness(c) %Eq. 25
            v2=x(b,:)-x(c,:);
        else
            v2=x(c,:)-x(b,:);
        end
        
        if rand<=rand
            
            u=1/3*(x(i,:)+bestSol+mo); %Eq . 19
            deta=sqrt(1/3*((x(i,:)-u).^2 ...
                    +(bestSol-u).^2+(mo-u).^2)); %Eq. 20
            
            vc1=rand(1,d);
            vc2=rand(1,d);
            
            %Eq. 21 ////////////////////////////////////////////
            Z1=sqrt(-1*log(vc2)).*cos(2*pi.*vc1);
            Z2=sqrt(-1*log(vc2)).*cos((2*pi.*vc1)+pi);
            a = rand;
            b = rand;
            if a<=b
                eta = (u+deta.*Z1);
            else
                eta = (u+deta.*Z2);
            end
            
            newsol = eta;
            %///////////////////////////////////////////////////////
        else
            beta=rand;
            v = x(i,:) +beta*abs(randn).*v1 ...   %Eq. 23
                                    +(1-beta)*abs(randn).*v2;
                                
            newsol = v;
        end
        
        newsol = max(newsol, lb);
        newsol = min(newsol, ub);
        
        newfitness =  obj(newsol);
        
        if newfitness<fitness(i)  %Eq. 27
            x(i,:) = newsol;
            if fitness(i) < bestFitness
                bestSol = x(i,:);
                bestFitness = fitness(i);
            end
        end
    end
    
    cgcurve(it)=bestFitness;
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(bestFitness,15)]);
end
 

end



