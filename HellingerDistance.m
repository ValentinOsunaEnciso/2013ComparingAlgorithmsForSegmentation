%% It calculates Hellinger distance, based on: Distance, by Herve Abdi %%%%
%%%   Valentin Osuna-Enciso, CIC-IPN, Abril, 2012 %%%%%%%%%%%%%%%%%%%%%%%%%
function d=HellingerDistance(A,B)
rA=sqrt(A);
rB=sqrt(B);
d=(rA-rB).^2;
d=sum(d);
d=sqrt(d);
