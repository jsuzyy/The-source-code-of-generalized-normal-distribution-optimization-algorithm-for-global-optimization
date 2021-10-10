%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA111
% Project Title: Implementation of TLBO in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
%
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
%
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%


function [BestCost,BestValue]=TLBO(fhd,nPop,nVar,VarMin,VarMax,MaxIt,X)
%% Problem Definition
MaxIt=MaxIt/(2*nPop);
% Cost Function


%nVar = 10;          % Number of Unknown Variables
VarSize = [1 nVar]; % Unknown Variables Matrix Size

%VarMin = -10;       % Unknown Variables Lower Bound
%VarMax =  10;       % Unknown Variables Upper Bound

%% TLBO Parameters

%MaxIt = 1000;        % Maximum Number of Iterations

%nPop = 50;           % Population Size
%% Initialization
% Empty Structure for Individuals

% Initialize Best Solution
BestSol = inf;
% Initialize Population Members
for i=1:nPop
    pop(i,:) = X(i,:);
    pop_Cost(i) = fhd(X(i,:));
    
    if  pop_Cost(i)  < BestSol
        BestSol =  pop_Cost(i) ;
        Best=pop(i,:);
    end
end

% Initialize Best Cost Record
BestCosts = zeros(MaxIt,1);
BestCost(1)=BestSol ;
%% TLBO Main Loop

for it=2:MaxIt
    
    % Calculate Population Mean
    Mean = 0;
    for i=1:nPop
        Mean = Mean + pop(i,:);
    end
    Mean = Mean/nPop;
    [~,la]=min(pop_Cost);
    [~,lb]=max(pop_Cost);
    % Select Teacher
    Teacher = pop(la,:);
    worst=pop(lb,:);
    
    % Teacher Phase
    for i=1:nPop
        % Create Empty Solution
        
        
        % Teaching Factor
        TF = randi([1 2]);
        
        % Teaching (moving towards teacher)
        newsol = pop(i,:)+ rand(VarSize).*(Teacher - TF*Mean);
        % Clipping
        newsol = max(newsol, VarMin);
        newsol= min(newsol, VarMax);
        
        % Evaluation
        newsol_Cost = fhd(newsol);
        
        % Comparision
        if newsol_Cost<pop_Cost(i)
            pop_Cost(i) = newsol_Cost;
            pop(i,:)=newsol;
            if pop_Cost(i) < BestSol
                BestSol = pop_Cost(i);
                Best=pop(i,:);
            end
        end
    end
    
    % Learner Phase
    for i=1:nPop
        
        A = 1:nPop;
        A(i)=[];
        j = A(randi(nPop-1));
        
        Step = pop(i,:) - pop(j,:);
        if pop_Cost(j) < pop_Cost(i)
            Step = -Step;
        end
        
        % Create Empty Solution
        
        % Teaching (moving towards teacher)
        
        newsol = pop(i,:) + rand(VarSize).*Step;
        
        % Clipping
        newsol = max(newsol, VarMin);
        newsol = min(newsol, VarMax);
        
        % Evaluation
        newsol_Cost = fhd(newsol);
        
        % Comparision
        if newsol_Cost<pop_Cost(i)
            pop_Cost(i) = newsol_Cost;
            pop(i,:)=newsol;
            if pop_Cost(i) < BestSol
                BestSol = pop_Cost(i);
                Best=pop(i,:);
            end
        end
    end
    
    % Store Record for Current Iteration
    BestCosts(it) = BestSol;
    
    % Show Iteration Information
    %  disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    BestCost(it)=BestSol;
    
    %   disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestSol.Cost,15)]);
end
BestValue=BestSol;
%% Results

end
