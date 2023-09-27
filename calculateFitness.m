function fFitness=calculateFitness(F_x_)
fFitness=zeros(size(F_x_));
ind=find(F_x_>=0);
fFitness(ind)=1./(F_x_(ind)+1);
ind=find(F_x_<0);
fFitness(ind)=1+abs(F_x_(ind));
