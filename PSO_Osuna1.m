%% PSO simple+Thresholding,Método Suma de Gaussianas, Osuna, CIC, 2011
function [temp k EvaluacioneFO tiempo]=PSO_MGF
clear all
format long
DB=imread('Im003_1.jpg'); 
numClases=3;    %Numero de clases que deseo obtener(a-priori)
DB=rgb2gray(DB);%Convierto a escala de grises
%DB=DB(:,:,1);
H=imhist(DB);   %Calculo del Histograma
H=H/sum(H);     %Se normaliza Histograma experimental(suma de Hi=1)
Amax=max(H);    %Alturas maximas
%% INICIALIZACION PSO: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmax=2000; k=1;  %Nmax=Numero maximo de iteraciones, k=Iteracion actual
particulas=90;  % Numero de particulas del enjambre
dimensiones=numClases*3;  % Num de parametros de funciones Gaussianas(alt,med,desv)
L=size(H,1);    % Numero de niveles de gris 
rho=4.6;wt=3; c1=2;c2=2;r1=0.48;r2=0.48; chi=0.72968; %rho=parametroDeControl, wo=Peso inicial
x_high=[Amax;L-1;(L-1)/12;Amax;L-1;(L-1)/12;Amax;L-1;(L-1)/12];
x_low=[0;0;0;0;0;0;0;0;0];%[Amax/10;1;1;Amax/10;1;1;Amax/10;1;1;];
v_high=[Amax;L-1;(L-1)/12;Amax;L-1;(L-1)/12;Amax;L-1;(L-1)/12];
v_low=[0;0;0;0;0;0;0;0;0];
for ind1=1:particulas
    for ind2=1:dimensiones 
        V(ind1,ind2)=v_low(ind2,1)+rand()*(v_high(ind2,1)-v_low(ind2,1));
        X(ind1,ind2)=x_low(ind2,1)+rand()*(x_high(ind2,1)-x_low(ind2,1)); 
    end
end
atorado=0;
pbest=X;gbest=X(1,:);f_pbest=99*ones(particulas,1);f_gbest=99;
% rng('shuffle');
tic
%% 2. REPETIR: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while k<=Nmax && f_gbest(1)>0.1186 && atorado<200
    %% 3. EVALUAR CADA PARTICULA,COMPARAR APTITUDES GLOBAL Y POR PARTICULA:
    tempo=f_gbest(1);
    [f_pbest,pbest,f_gbest,gbest]=sumgauss(X,H,L,particulas,dimensiones,...
                                    f_pbest,pbest,f_gbest,gbest); 
    if(tempo==f_gbest(1))
        atorado=atorado+1;
    else
        atorado=0;
    end
    %% 4. CAMBIAR LA POSICION Y VELOCIDAD DE CADA PARTICULA: %%%%%%%%%%%%%%
    for ind1=1:particulas
        %%Version:MultilevelThresholdingAlgorithmBasedOnPSOforImageSegmenta
        %r1=4.0*r1*(1-r1);
        %r2=4.0*r2*(1-r2);
        %wt=wt*exp(-rho*k/Nmax);
        %%Version Original PSO
        %r1=rand()*c1; r2=rand()*c2;
        wt=randn()*k/Nmax;
        for ind2=1:dimensiones
            r1=rand()*c1; r2=rand()*c2;
            V(ind1,ind2)=wt*V(ind1,ind2)+...
            r1*(pbest(ind1,ind2)-X(ind1,ind2))+...
            r2*(gbest(1,ind2)-X(ind1,ind2));
            %if V(ind1,ind2)>v_high(ind2,1)% ||V(ind1,ind2)<v_low(ind2,1)
            %    V(ind1,ind2)=v_high(ind2,1);
            %end
            %if V(ind1,ind2)<v_low(ind2,1)
            %    V(ind1,ind2)=v_low(ind2,1);
            %end
            X(ind1,ind2)=X(ind1,ind2)+V(ind1,ind2); %Posicion
            if X(ind1,ind2)>x_high(ind2,1) %X(ind1,ind2)<x_low(ind2,1)
                X(ind1,ind2)=x_high(ind2,1); 
            end
            if X(ind1,ind2)<x_low(ind2,1)
                X(ind1,ind2)=x_low(ind2,1); 
            end
        end
    end
    disp(sprintf('Iterac=%d, Fitnes=%f\n',k,f_gbest(1)));
    k=k+1;
end
tiempo=toc;
temp=f_gbest(1);
grafica(gbest, dimensiones,H,DB);
end
%%Termina PSO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EVALUAR FUNCION OTSU, COMPARAR APTITUDES PERSONALES Y GLOBALES: %%%%%%%%
function [f_pbest,pbest,f_gbest,gbest]=otsu(X,H,L,P,D,f_pbest,pbest,f_gbest,gbest)
  if D==1
    for ind2=1:P
        t1=X(ind2,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        q1=sum(H(1:t1,1));
        if q1==0
            temp1=0;
        else
            m1=round(sum((1:t1)*H(1:t1,1))/q1);           
            var1=sum((((1:t1)-m1).^2)*H(1:t1,1)/q1);
            q2=sum(H(t1+1:L,1));
            if q2==0
                temp1=0;
            else
                m2=round(sum((t1+1:L)*H(t1+1:L,1))/q2);        
                var2=sum((((t1+1:L)-m2).^2)*H(t1+1:L,1)/q2);      
                mt=q1*m1+q2*m2;
                temp1=(q1*(m1-mt)^2)+(q2*(m2-mt)^2);%Maximizar Var_bc
                %temp1=q1*var1+q2*var2;             %Minimizar Var_wc                 
            end   
        end
        if temp1>f_pbest(ind2,1)
            f_pbest(ind2,1)=temp1;
            pbest(ind2,1)=X(ind2,:);
            if temp1>f_gbest
                f_gbest=temp1;
                gbest=X(ind2,:);
            end
        end
    end
  else
    for ind2=1:P
        t1=X(ind2,1);
        t2=X(ind2,2);
        if(t1>t2)
            t1=X(ind2,2);t2=X(ind2,1);
            X(ind2,2)=t2;X(ind2,1)=t1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        q1=sum(H(1:t1,1));
        if q1==0
            temp1=0;
        else
            m1=round(sum((1:t1)*H(1:t1,1))/q1);           
            var1=sum((((1:t1)-m1).^2)*H(1:t1,1)/q1);
            q2=sum(H(t1+1:t2,1));
            if q2==0
                temp1=0;
            else
                m2=round(sum((t1+1:t2)*H(t1+1:t2,1))/q2);        
                var2=sum((((t1+1:t2)-m2).^2)*H(t1+1:t2,1)/q2);
                q3=sum(H(t2+1:L,1));
                if q3==0
                    temp1=0;
                else
                    m3=round(sum((t2+1:L)*H(t2+1:L,1))/q3);        
                    var3=sum((((t2+1:L)-m2).^2)*H(t2+1:L,1)/q3);
                    mt=q1*m1+q2*m2+q3*m3;
                    temp1=(q1*(m1-mt)^2)+(q2*(m2-mt)^2)+(q3*(m3-mt)^2);%Maximizar Var_bc
                    %temp1=q1*var1+q2*var2;             %Minimizar Var_wc
                end
            end   
        end
        if temp1>f_pbest(ind2,1)
            f_pbest(ind2,1)=temp1;
            pbest(ind2,:)=X(ind2,:);
            if temp1>f_gbest
                f_gbest=temp1;
                gbest=X(ind2,:);
            end
        end
    end    
  end
end
%%
function error=MGF(x,H)
    xp=0:1:255;
    mix(:,1)=(x(1,1)*exp(-((xp-round(x(4,1))).^2)...
   /(2*(x(7,1)^2))))+(x(2,1)*exp(-((xp-round(x(5,1))).^2)...
   /(2*(x(8,1)^2))))+(x(3,1)*exp(-((xp-round(x(6,1))).^2)...
   /(2*(x(9,1)^2))));
%     medias=[round(x(4,1)),round(x(5,1)),round(x(4,1))];
%     medias=sort(medias);
%     error3=HellingerDistance(H,mix)+(3/((medias(3)-medias(2))+(medias(2)-medias(1))+eps));
    error=HellingerDistance(H,mix);
end
%% EVALUAR FUNCION GAUSSIANAS, COMPARAR APTITUDES PERSONALES Y GLOBALES: %%
function [f_pbest,pbest,f_gbest,gbest]=sumgauss(X,H,L,P,D,f_pbest,pbest,f_gbest,gbest)
%EVALUACION DE LA FUNCION DE COSTO: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xp=0:1:255;
  if D==9    
    for ind1=1:P
        gaussianas(:,ind1)=...
        (X(ind1,1)*exp(-((xp-round(X(ind1,2))).^2)/(2*(X(ind1,3)^2))))+...
        (X(ind1,4)*exp(-((xp-round(X(ind1,5))).^2)/(2*(X(ind1,6)^2))))+...
        (X(ind1,7)*exp(-((xp-round(X(ind1,8))).^2)/(2*(X(ind1,9)^2))));                           
        temp1=HellingerDistance(H,gaussianas(:,ind1));
        if temp1<f_pbest(ind1,1)
            f_pbest(ind1,1)=temp1;
            pbest(ind1,:)=X(ind1,:);
            if temp1<f_gbest               
                f_gbest=temp1;
                gbest=X(ind1,:);
            end
        end
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
%Recibe gbest: La mejor partícula, D: dimensiones de cada particula
%H: histograma de la imagen, DB: imagen en escala de gris
function grafica(gbest,dimensiones,H,DB)
  xp=0:1:255;
  if dimensiones==9    
    valM1=round(gbest(1,2));valM2=round(gbest(1,5));valM3=round(gbest(1,8));
    valA1=gbest(1,1);valA2=gbest(1,4);valA3=gbest(1,7);
    valDE1=(gbest(1,3));valDE2=(gbest(1,6));valDE3=(gbest(1,9));
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
    Resultado=(A1*exp(-((xp-M1).^2)/(2*(DE1^2))))+(A2*exp(-((xp-M2).^2)/(2*(DE2^2))))+(A3*exp(-((xp-M3).^2)/(2*(DE3^2))));
    plot(Resultado,'k')%,Hold on,plot((ampli1*exp(-((x-media1).^2)/(2*(DE1^2))))),plot((ampli2*exp(-((x-media2).^2)/(2*(DE2^2))))),plot((ampli3*exp(-((x-media3).^2)/(2*(DE3^2)))))
    hold on
    plot(H,'k','LineWidth',2)
    figure
    plot((A3*exp(-((xp-M3).^2)/(2*(DE3^2)))),'k--','LineWidth',2),hold on
    plot((A2*exp(-((xp-M2).^2)/(2*(DE2^2)))),'k-.','LineWidth',2)
    plot((A1*exp(-((xp-M1).^2)/(2*(DE1^2)))),'k','LineWidth',2)
    %plot(Resultado,'g','LineWidth',2),title('Resultado')
%%%%%Realizo umbralizacion imagen escala de grises:%%%%%%%%%%%%%%%%%%%%%%%
    a1=(DE1^2)-(DE2^2);
    a2=(DE2^2)-(DE3^2);
    b1=2*((M1*(DE2^2))-(M2*(DE1^2)));
    b2=2*((M2*(DE3^2))-(M3*(DE2^2)));
    c1=((DE1*M2)^2)-((DE2*M1)^2)+(2*((DE1*DE2)^2)*log((DE2*A1)/(DE1*A2)));
    c2=((DE2*M3)^2)-((DE3*M2)^2)+(2*((DE3*DE2)^2)*log((DE3*A2)/(DE2*A3)));
    T1a=(-b1+sqrt((b1^2)-(4*a1*c1)))/(2*a1);
    T1b=(-b1-sqrt((b1^2)-(4*a1*c1)))/(2*a1);
    T2a=(-b2+sqrt((b2^2)-(4*a2*c2)))/(2*a2);
    T2b=(-b2-sqrt((b2^2)-(4*a2*c2)))/(2*a2);
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
  elseif dimensiones==6
    valM1=round(gbest(1,2));valM2=round(gbest(1,5));
    valA1=gbest(1,1);valA2=gbest(1,4);
    valDE1=(gbest(1,3));valDE2=(gbest(1,6));
    [M ind1]=sort([valM1,valM2]);
    if(ind1(1)==1)
        DE1=valDE1;A1=valA1;M1=valM1;
    else
        DE1=valDE2;A1=valA2;M1=valM2;
    end
    if(ind1(2)==1)
        DE2=valDE1;A2=valA1;M2=valM1;
    else
        DE2=valDE2;A2=valA2;M2=valM2;
    end
    Resultado=(A1*exp(-((xp-M1).^2)/(2*(DE1^2))))+...
        (A2*exp(-((xp-M2).^2)/(2*(DE2^2))));
% plot(Resultado,'k')%,Hold on,plot((ampli1*exp(-((x-media1).^2)/(2*(DE1^2))))),plot((ampli2*exp(-((x-media2).^2)/(2*(DE2^2))))),plot((ampli3*exp(-((x-media3).^2)/(2*(DE3^2)))))
% Hold on
    plot(H,'r','LineWidth',2),hold on
%figure
    plot((A2*exp(-((xp-M2).^2)/(2*(DE2^2)))),'k--','LineWidth',2),hold on
    plot((A1*exp(-((xp-M1).^2)/(2*(DE1^2)))),'k','LineWidth',2)
    plot(Resultado,'g','LineWidth',2),title('Resultado')
%%%%%Realizo umbralizacion imagen escala de grises:%%%%%%%%%%%%%%%%%%%%%%%
    a1=(DE1^2)-(DE2^2);
    b1=2*((M1*(DE2^2))-(M2*(DE1^2)));
    c1=((DE1*M2)^2)-((DE2*M1)^2)+(2*((DE1*DE2)^2)*log((DE2*A1)/(DE1*A2)));
    T1a=(-b1+sqrt((b1^2)-(4*a1*c1)))/(2*a1);
    T1b=(-b1-sqrt((b1^2)-(4*a1*c1)))/(2*a1);
    [fila columna]=size(DB);
    for ind1=1:fila
        for ind2=1:columna
            if (DB(ind1,ind2)<=T1b)&&(DB(ind1,ind2)>=0)
                DBsegmented(ind1,ind2)=0;
            else
                DBsegmented(ind1,ind2)=1;
            end
        end
    end
    figure,DBsegmented=mat2gray(DBsegmented);
    imsHow(DBsegmented)
  end
end