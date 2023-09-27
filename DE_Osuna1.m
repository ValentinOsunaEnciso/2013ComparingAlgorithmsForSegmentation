%DE1 simple+Thresholding, Osuna, CIC, 2011
%A novel multithreshold segmentation approach based on DE optimization,'09
function temp=DE_Osuna1
clear all
format long
DB=imread('Im003_1.jpg'); 
numClases=3;    %Numero de clases que deseo obtener(a-priori)
DB=rgb2gray(DB);%Convierto a escala de grises
H=imhist(DB);   %Calculo del Histograma
H=H/sum(H); Amax=max(H);%Se normaliza Histograma experimental(suma de Hi=1)
L=size(H,1);
w=1.5;          %Factor de penalizaci�n
x_high=[Amax;Amax;Amax;L-1;L-1;L-1;(L-1)/12;(L-1)/12;(L-1)/12];
x_low=[0;0;0;0;0;0;0;0;0];x=[]; xp=0:1:255;
%INICIALIZACION DE1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np=90;          %Tamano de la poblacion.90
F=0.25;         %Factor de escalamiento
Cr=0.8;         %Probabilidad de cruza
D=numClases*3;  %Para 3 clases, 9 dimensiones.
gmax=2000;       %Numero maximo de iteraciones
t=0;            %Contador de iteraciones
for ind1=1:Np
    for ind2=1:D
       x(ind2,ind1)=x_low(ind2,1)+rand()*(x_high(ind2,1)-x_low(ind2,1));  
    end
end
%temp=x(:,ind(1:end));x=temp; %Ordeno en base a los fitness
f_x(1,1)=1;
% rng('shuffle');
while (t<gmax && f_x(1,1)>0.0000021)
  if D==9
   Ampli1=x(1,1:Np); Ampli2=x(2,1:Np); Ampli3=x(3,1:Np);
   Med1=round(x(4,1:Np)); Med2 = round(x(5,1:Np)); Med3 = round(x(6,1:Np));
   desvest1 =x(7,1:Np);  desvest2 = x(8,1:Np); desvest3 = x(9,1:Np);
   %EVALUACIÓN DE LA FUNCIÓN DE COSTO: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for ind1=1:Np
   gaussianas(:,ind1)=(Ampli1(1,ind1)*exp(-((xp-Med1(1,ind1)).^2)...
   /(2*(desvest1(1,ind1)^2))))+(Ampli2(1,ind1)*exp(-((xp-Med2(1,ind1)).^2)...
   /(2*(desvest2(1,ind1)^2))))+(Ampli3(1,ind1)*exp(-((xp-Med3(1,ind1)).^2)...
   /(2*(desvest3(1,ind1)^2))));
      errortotal(:,ind1)=((H-gaussianas(:,ind1))).^2;
      error(ind1)=sum(errortotal(:,ind1))/L;
   end
   [f_x ind]=sort(error);%Ordeno deMenor aMayor la evaluacion de la funcion
   x_best=x(:,ind(1));%Guardo el mejor;Van Np Evaluaciones De La Funcion.
   for ind1=1:Np
       r1=randi(Np);r2=randi(Np);
       while(r1==r2==ind1)
          r1=randi(Np); r2=randi(Np);%Generados sean diferentes.
       end
       v(:,1)=x_best+F*(x(:,r1)-x(:,r2));%Mutacion
       u(:,1)=x(:,ind1);
       j_rand=randi(D);
       for ind2=1:D                         %Genero vector de prueba
          if(rand()<Cr || ind2==j_rand)     %Cruza
             u(j_rand,1)=v(j_rand,1);
          end
       end
       temp1(:,1)=(u(1,1)*exp(-((xp-round(u(4,1))).^2)...
   /(2*(u(7,1)^2))))+(u(2,1)*exp(-((xp-round(u(5,1))).^2)...
   /(2*(u(8,1)^2))))+(u(3,1)*exp(-((xp-round(u(6,1))).^2)...
   /(2*(u(9,1)^2))));
      errortotal(:,1)=((H-temp1)).^2;
      temp=sum(errortotal(:,ind1))/L;
      if(temp<f_x(ind(ind1)))
          x(:,ind(ind1))=u(:,1);
      end
   end
   x(:,ind(1))=x_best;
  end
t=t+1;
disp(sprintf('Iterac=%d, Fitnes=%f, a1=%f, m1=%f, d1=%f\n',t,f_x(1,1),x_best(1,1),x_best(2,1),x_best(3,1)));
end
grafica(x_best,D,H,DB)
temp=f_x(1,1);
end
%%Termina DE1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCION GRAFICA GAUSSIANAS E IMAGEN SEGMENTADA: %%%%%%%%%%%%%%%%%%%%%%%%
%Recibe gbest: La mejor partícula, D: dimensiones de cada particula
%H: histograma de la imagen, DB: imagen en escala de gris
function grafica(x_best,D,H,DB)
  if D==9
    xp=0:1:255;
    valM1=round(x_best(4,1));valM2=round(x_best(5,1));valM3=round(x_best(6,1));
    valA1=x_best(1,1);valA2=x_best(2,1);valA3=x_best(3,1);
    valDE1=x_best(7,1);valDE2=x_best(8,1);valDE3=x_best(9,1);
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
    plot(Resultado,'k--','LineWidth',2)%,Hold on,plot((ampli1*exp(-((x-media1).^2)/(2*(DE1^2))))),plot((ampli2*exp(-((x-media2).^2)/(2*(DE2^2))))),plot((ampli3*exp(-((x-media3).^2)/(2*(DE3^2)))))
    hold on
    plot(H,'r','LineWidth',2),figure
    plot((A3*exp(-((xp-M3).^2)/(2*(DE3^2)))),'k--','LineWidth',2),hold on
    plot((A2*exp(-((xp-M2).^2)/(2*(DE2^2)))),'k-.','LineWidth',2)
    plot((A1*exp(-((xp-M1).^2)/(2*(DE1^2)))),'k','LineWidth',2)
    %plot(Resultado,'k'),title('Resultado')
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
            elseif (DB(ind1,ind2,1)<=T2b)&&(DB(ind1,ind2)>T1b)
                DBsegmented(ind1,ind2)=0.5;
            elseif(DB(ind1,ind2,1)>T2b)
                DBsegmented(ind1,ind2)=1;
            end
        end
    end
    figure,DBsegmented=mat2gray(DBsegmented);
    imshow(DBsegmented)
  end
end
