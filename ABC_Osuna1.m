function temp=ABC_Osuna1
%/* ABC algorithm coded using MATLAB language */
%/* Artificial Bee Colony (ABC) is one of the most recently defined algorit
%hms by Dervis Karaboga in 2005, motivated by the intelligent behavior of 
%honey bees. Referance Papers:
%1)D. Karaboga,AN IDEA BASED ON HONEY BEE SWARM FOR NUMERICAL OPTIMIZATION,
%TECHNICAL REPORT-TR06, Erciyes University, Engineering Faculty, Computer 
%Engineering Department 2005. 2)D. Karaboga, B. Basturk, A powerful and 
%Efficient Algorithm for Numerical Function Optimization: Artificial Bee 
%Colony (ABC) Algorithm, Journal of Global Optimization, Volume:39,Issue:3,
%pp:459-171, November 2007,ISSN:0925-5001 , doi: 10.1007/s10898-007-9149-x 
%3)D. Karaboga, B. Basturk, On The Performance Of Artificial Bee Colony 
%(ABC) Algorithm, Applied Soft Computing,Volume 8, Issue 1, January 2008, 
%Pages 687-697. 4)D. Karaboga, B. Akay, A Comparative Study of Artificial 
%Bee Colony Algorithm,  Applied Mathematics and Computation, 214, 108-132, 
%2009.
%Copyright,2009,ErciyesUniversity,IntelligentSystemsResearchGroup,TheDept. 
%of Computer Engineering. Contact:Dervis Karaboga (karaboga@erciyes.edu.tr)
%Bahriye Basturk Akay (bahriye@erciyes.edu.tr)
clear all
close all
clc
DB=imread('Im003_1.jpg'); 
numClases=3;    %Numero de clases que deseo obtener(a-priori)
DB=rgb2gray(DB);%Convierto a escala de grises
%DB=DB(:,:,1);
H=imhist(DB);   %Calculo del Histograma
H=H/sum(H);     %Se normaliza Histograma experimental(suma de Hi=1)
Amax=max(H);    %Alturas m√°ximas
L=size(H,1);    % Numero de niveles de gris 
% CONTROL PARAMETERS OF ABC ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NP=90;          %The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2;%TheNumberOfFoodSourcesEqualsTheHalf of the colony size
%AFoodSourceWhichCouldn'tBeImprovedThrough"limit"trialsIsAbandonedByIts
limit=10; %employed bee.
Nmax=2000; %/*The number of cycles for foraging {a stopping criteria}*/
%/* Problem specific variables*/
%objfun='Sphere'; %cost function to be optimized
D=numClases*3;  % Num de parametros de funciones Gaussianas(alt,med,desv)
x_high=[Amax;L-1;(L-1)/12;Amax;L-1;(L-1)/12;Amax;L-1;(L-1)/12];
x_low=[0;0;0;0;0;0;0;0;0];%[Amax/10;1;1;Amax/10;1;1;Amax/10;1;1;];
for ind1=1:FoodNumber
    for ind2=1:D 
        Foods(ind1,ind2)=(x_low(ind2,1)+rand()*...
           (x_high(ind2,1)-x_low(ind2,1))); 
    end
end
%Foods [FoodNumber][D]; /*Foods is the population of food sources. 
%Each row of Foods matrix is a vector holding D parameters to be optimized. 
%The number of rows of Foods matrix equals to the FoodNumber
%F_x_[FoodNumber];  /*f is a vector holding objective function values 
%associated with food sources. Fitness[FoodNumber]; fitness is a vector 
%holding fitness (quality) values associated with food sources*/
%trial[FoodNumber]; /*trial is a vector holding trial numbers through which
%solutions can not be improved. prob[FoodNumber]; /*prob is a vector 
%holding probabilities of food sources (solutions) to be chosen*/
%solution [D]; /*New solution (neighbour) produced by 
%v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter 
%and k is a randomlu chosen solution different from i.ObjValSol; Objective 
%function value of new solution. FitnessSol; Fitness value of new solution
%neighbour, param2change; /*param2change corrresponds to j, neighbour 
%corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
%GlobalMin; /*Optimum solution obtained by ABC algorithm*/
%GlobalParams[D]; /*Parameters of the optimum solution*/
%GlobalMins[runtime];GlobalMins holdsThe GlobalMin ofEach runInMultipleRuns
%GlobalMins=zeros(1,runtime); 
% /*All food sources are initialized */
%/*Variables are initialized in the range [lb,ub]. If each parameter has 
%different range, use arrays lb[j], ub[j] instead of lb and ub 
% Range = repmat((ub-lb),[FoodNumber 1]);
% Lower = repmat(lb, [FoodNumber 1]);
% Foods = rand(FoodNumber,D) .* Range + Lower;
%F_x_=feval(objfun,Foods);
F_x_=sumgauss(Foods,H,L,FoodNumber,D);
Fitness=calculateFitness(F_x_*10000);
%reset trial counters
trial=zeros(1,FoodNumber);
%/*The best food source is memorized*/
BestInd=find(F_x_==min(F_x_));
BestInd=BestInd(end);
GlobalMin=F_x_(BestInd);
GlobalParams=Foods(BestInd,:);
k=1;
% rng('shuffle');
while (k <= Nmax && GlobalMin>0.000001)
    %%%%%%%%% EMPLOYED BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%
    for ind1=1:(FoodNumber)        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;        
       %ArandomlyChosenSolutionIsUsedInProducingaMutantSolutionOfsolution i
       neighbour=fix(rand*(FoodNumber))+1;       
        %/*Randomly selected solution must be different from the solution i       
            while(neighbour==ind1)
                neighbour=fix(rand*(FoodNumber))+1;
            end;        
       sol=Foods(ind1,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(ind1,Param2Change)+...
           (Foods(ind1,Param2Change)-Foods(neighbour,Param2Change))*...
           (rand-0.5)*2;        
       %IfGeneratedParameterValueIsOutOfBoundaries,ItIsShiftedOntoTheBounda
       ind=find(sol<x_low(:,1)');
       sol(ind)=x_low(ind,1)';
       ind=find(sol>x_high(:,1)');
       sol(ind)=x_high(ind,1)';        
       %Evaluate new solution
       ObjValSol=sumgauss(sol,H,L,1,D);%feval(objfun,sol);       
       FitnessSol=calculateFitness(ObjValSol*10000);        
       %GreedySelectionIsAppliedBetweenTheCurrentSolution i and its mutant
       %IfTheMutantSolutionIsBetterThanTheCurrent solution i, replace the 
       %solution with the mutant and reset the trial counter of solution i
       if (FitnessSol>Fitness(ind1)) %¬ø>, o <? Originalmente >
          Foods(ind1,:)=sol;
          Fitness(ind1)=FitnessSol;
          F_x_(ind1)=ObjValSol;
          trial(ind1)=0;
       else
          %Ifsolution i can'tBeImproved,IncreaseTrialCounter  
          trial(ind1)=trial(ind1)+1;
       end;                  
    end;
    % CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A food source is chosen with the probability which is proportioal to 
    % its quality Different schemes can be used to calculate the 
    % probability values. For example prob(i)=fitness(i)/sum(fitness) or
    % in a way used in the metod below prob(i)=a*fitness(i)/max(fitness)+b
    % probability values are calculated by using fitness values and 
    % normalized by dividing maximum fitness value*/
    prob=(0.9.*Fitness./max(Fitness))+0.1;  
    % ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind1=1;
    t=0;
    while(t<FoodNumber)
        if(rand<prob(ind1))
            t=t+1;
            %/*The parameter to be changed is determined randomly*/
            Param2Change=fix(rand*D)+1;        
    %ArandomlyChosenSolutionIsUsedInProducingaMutantSolutionOfTheSolution i
            neighbour=fix(rand*(FoodNumber))+1;       
            %RandomlySelectedSolutionMustBeDifferentFromTheSolution i        
            while(neighbour==ind1)
                neighbour=fix(rand*(FoodNumber))+1;
            end;        
            sol=Foods(ind1,:);
            %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
            sol(Param2Change)=Foods(ind1,Param2Change)+...
                (Foods(ind1,Param2Change)-Foods(neighbour,Param2Change))*...
                (rand-0.5)*2;        
       %IfGeneratedParameterValueIsOutOfBoundaries,ItIsShiftedOntoTheBounda
        ind=find(sol<x_low(:,1)');
        sol(ind)=x_low(ind,1)';
        ind=find(sol>x_high(:,1)');
        sol(ind)=x_high(ind,1)';        
        %Evaluate new solution
        ObjValSol=sumgauss(sol,H,L,1,D);%feval(objfun,sol);       
        FitnessSol=calculateFitness(ObjValSol*10000);        
        %GreedySelectionIsAppliedBetweenTheCurrentSolution i and its mutant
        %If the mutant solution is better than the current solution i, 
        %replaceTheSolutionWithTheMutantAndResetTheTrialCounterOfSolution i
            if (FitnessSol>Fitness(ind1))
                Foods(ind1,:)=sol;
                Fitness(ind1)=FitnessSol;
                F_x_(ind1)=ObjValSol;
                trial(ind1)=0;
            else
                %IfSolution ind1 can'tBeImproved,IncreaseIts trial counter
                trial(ind1)=trial(ind1)+1; 
            end;
        end;    
        ind1=ind1+1;
        if (ind1==(FoodNumber)+1) 
            ind1=1;
        end;   
    end; 
    %/*The best food source is memorized*/
    ind=find(F_x_==min(F_x_));
    ind=ind(end);
    if (F_x_(ind)<GlobalMin)
         GlobalMin=F_x_(ind);
         GlobalParams=Foods(ind,:);
    end;       
    %%%%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DetermineThe food sources whose trial counter exceeds the "limit"value 
    %In Basic ABC, only one scout is allowed to occur in each cycle
    ind=find(trial==max(trial));
    ind=ind(end);
    if (trial(ind)>limit)
        trial(ind)=0;
        %sol=(ub-lb).*rand(1,D)+lb;
        for ind2=1:D 
            sol(1,ind2)=(x_low(ind2,1)+rand()*...
           (x_high(ind2,1)-x_low(ind2,1))); 
        end
        ObjValSol=sumgauss(sol,H,L,1,D);%feval(objfun,sol);       
        FitnessSol=calculateFitness(ObjValSol*10000); 
        Foods(ind,:)=sol;
        Fitness(ind)=FitnessSol;
        F_x_(ind)=ObjValSol;
    end;
   fprintf('Iter=%d F_x_=%g\n',k,GlobalMin);
    k=k+1;
end % End of ABC
%GlobalMins(r)=GlobalMin;
grafica(GlobalParams', D,H,DB);
temp=GlobalMin;
end

%% EVALUAR FUNCION GAUSSIANAS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp2=sumgauss(X,H,L,P,D)
%EVALUACI”N DE LA FUNCI”N DE COSTO: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xp=0:1:255;
  if D==9    
    for ind1=1:P
        gaussianas(:,ind1)=...
        (X(ind1,1)*exp(-((xp-round(X(ind1,2))).^2)/(2*(X(ind1,3)^2))))+...
        (X(ind1,4)*exp(-((xp-round(X(ind1,5))).^2)/(2*(X(ind1,6)^2))))+...
        (X(ind1,7)*exp(-((xp-round(X(ind1,8))).^2)/(2*(X(ind1,9)^2))));            
        temp1=((H-gaussianas(:,ind1))).^2;
        %temp2=.5*abs(sum(gaussianas(:,ind1))-1);
        temp2(:,ind1)=(sum(temp1))/L;        
    end
  elseif D==6
    for ind1=1:P
        gaussianas(:,ind1)=...
        (X(ind1,1)*exp(-((xp-round(X(ind1,2))).^2)/(2*(X(ind1,3)^2))))+...
        (X(ind1,4)*exp(-((xp-round(X(ind1,5))).^2)/(2*(X(ind1,6)^2))));            
        temp1=((H-gaussianas(:,ind1))).^2;
        %temp2=.5*abs(sum(gaussianas(:,ind1))-1);
        temp1=(sum(temp1))/L;        
        if temp1<f_pbest(ind1,1)
            f_pbest(ind1,1)=temp1;
            pbest(ind1,:)=X(ind1,:);
            if temp1<f_gbest               
                f_gbest=temp1;
                gbest=X(ind1,:);
            end
        end
    end  
  end
end

%% FUNCION GRAFICA GAUSSIANAS E IMAGEN SEGMENTADA: %%%%%%%%%%%%%%%%%%%%%%%%
%Recibe gbest: La mejor partÌcula, D: dimensiones de cada particula
%H: histograma de la imagen, DB: imagen en escala de gris
function grafica(x_best,D,H,DB)
  if D==9
    xp=0:1:255;
    valM1=round(x_best(2,1));valM2=round(x_best(5,1));valM3=round(x_best(8,1));
    valA1=x_best(1,1);valA2=x_best(4,1);valA3=x_best(7,1);
    valDE1=x_best(3,1);valDE2=x_best(6,1);valDE3=x_best(9,1);
    [M ind1]=sort([valM1,valM2,valM3]);
    if(ind1(1)==1)
        DE1=valDE1;A1=valA1;M1=valM1;
    elseif(ind1(1)==2)
        DE1=valDE2;A1=valA2;M1=valM2;
    else
        DE1=valDE3;A1=valA3;M1=valM3;
    end
    if(ind1(2)==1)
        DE2=valDE1;A2=valA1;M2=valM1;
    elseif(ind1(2)==2)
        DE2=valDE2;A2=valA2;M2=valM2;
    else
        DE2=valDE3;A2=valA3;M2=valM3;
    end
    if(ind1(3)==1)
        DE3=valDE1;A3=valA1;M3=valM1;
    elseif(ind1(3)==2)
        DE3=valDE2;A3=valA2;M3=valM2;
    else
        DE3=valDE3;A3=valA3;M3=valM3;
    end
    Resultado=(A1*exp(-((xp-M1).^2)/(2*(DE1^2))))+...
        (A2*exp(-((xp-M2).^2)/(2*(DE2^2))))+...
        (A3*exp(-((xp-M3).^2)/(2*(DE3^2))));
    plot(Resultado,'k')%,Hold on,plot((ampli1*exp(-((x-media1).^2)/(2*(DE1^2))))),plot((ampli2*exp(-((x-media2).^2)/(2*(DE2^2))))),plot((ampli3*exp(-((x-media3).^2)/(2*(DE3^2)))))
    hold on
    plot(H,'r'),figure
    plot((A3*exp(-((xp-M3).^2)/(2*(DE3^2)))),'k--'),hold on
    plot((A2*exp(-((xp-M2).^2)/(2*(DE2^2)))),'k-.')
    plot((A1*exp(-((xp-M1).^2)/(2*(DE1^2)))),'k')
    %plot(Resultado,'k'),title('Resultado')
%%%%%Realizo umbralizacion imagen escala de grises:%%%%%%%%%%%%%%%%%%%%%%%
    a1=(DE1^2)-(DE2^2);
    a2=(DE2^2)-(DE3^2);
    b1=2*((M1*(DE2^2))-(M2*(DE1^2)));
    b2=2*((M2*(DE3^2))-(M3*(DE2^2)));
    c1=((DE1*M2)^2)-((DE2*M1)^2)+(2*((DE1*DE2)^2)*log((DE2*A1)/(DE1*A2)));
    c2=((DE2*M3)^2)-((DE3*M2)^2)+(2*((DE3*DE2)^2)*log((DE3*A2)/(DE2*A3)));
    if a1~=0 &&a2~=0        
        T1a=(-b1+sqrt((b1^2)-(4*a1*c1)))/(2*a1);
        T1b=(-b1-sqrt((b1^2)-(4*a1*c1)))/(2*a1);
        T2a=(-b2+sqrt((b2^2)-(4*a2*c2)))/(2*a2);
        T2b=(-b2-sqrt((b2^2)-(4*a2*c2)))/(2*a2);
    else
        a1=1;a2=1;
        T1a=(-b1+sqrt((b1^2)-(4*a1*c1)))/(2*a1);
        T1b=(-b1-sqrt((b1^2)-(4*a1*c1)))/(2*a1);
        T2a=(-b2+sqrt((b2^2)-(4*a2*c2)))/(2*a2);
        T2b=(-b2-sqrt((b2^2)-(4*a2*c2)))/(2*a2);
    end
    [fila columna]=size(DB);
    for ind1=1:fila
        for ind2=1:columna
            if (DB(ind1,ind2)<=T1b)&&(DB(ind1,ind2)>=0)
                DBsegmented(ind1,ind2)=0;
            elseif (DB(ind1,ind2)<=T2b)&&(DB(ind1,ind2)>T1b)
                DBsegmented(ind1,ind2)=0.5;
            elseif(DB(ind1,ind2)>T2b)
                DBsegmented(ind1,ind2)=1;
            end
        end
    end
    figure,DBsegmented=mat2gray(DBsegmented);
    imshow(DBsegmented)
  end
end
