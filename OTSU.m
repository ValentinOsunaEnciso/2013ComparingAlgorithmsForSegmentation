% Otsu method for 2 thresholds, 3 classes; exhaustive version. %%%%%%%%%%%%
% Valentin Osuna-Enciso, CIC-IPN, Abril, 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It receives only the histogram being segmented. %%%%%%%%%%%%%%%%%%%%%%%%%
function [tiempo mhd]=OTSU
clear all
DB=imread('Im222_0.tif'); 
DB=rgb2gray(DB);%Convierto a escala de grises
H=imHist(DB);   %Calculo del Histograma
cont=1;
tic
for t1=1:255
    for t2=1:255
        q1=sum(H(1:t1,1))+eps;
        m1=floor(sum((0:t1-1)*H(1:t1,1))/q1); 
        d1=sum((((0:t1-1)-m1).^2)*H(1:t1,1)/q1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        q2=sum(H(t1+1:t2,1))+eps;
        m2=floor(sum((t1:t2-1)*H(t1+1:t2,1))/q2); 
        d2=sum((((t1:t2-1)-m2).^2)*H(t1+1:t2,1)/q2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        q3=sum(H(t2+1:end,1))+eps;
        m3=floor(sum((t2:255)*H(t2+1:end,1))/q3); 
        d3=sum((((t2:255)-m3).^2)*H(t2+1:end,1)/q3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp(cont)=q1*d1+q2*d2+q3*d3;   % Minimization
        umbrales(:,cont)=[t1;t2];
        cont=cont+1;
    end
end
tiempo=toc;
[~, temp]=min(temp);
umbra=umbrales(:,temp);T1b=umbra(1,1);T2b=umbra(2,1);
[fila columna]=size(DB);
    for ind1=1:fila
        for ind2=1:columna
            if (DB(ind1,ind2)<=T1b)&&(DB(ind1,ind2)>=0)
                DBsegmented(ind1,ind2)=0;
            elseif (DB(ind1,ind2,1)<=T2b)&&(DB(ind1,ind2)>T1b)
                DBsegmented(ind1,ind2)=1;
            elseif(DB(ind1,ind2,1)>T2b)
                DBsegmented(ind1,ind2)=1;
            end
        end
    end
    %figure,DBsegmented=mat2gray(DBsegmented);
    %imsHow(DBsegmented)
    %% Modified Hausdorff distance:
    AI=DBsegmented;
    BI=imread('Im222_0_GT.tif');BI=rgb2gray(BI(:,:,1:3));
    [AIbordes t]=edge(AI,'canny',.1);
    [BIbordes t]=edge(BI,'canny',.1);
    [A(:,1) A(:,2)]=find(AIbordes);
    [B(:,1) B(:,2)]=find(BIbordes);
    [ mhd ] = ModHausdorffDist( A, B );