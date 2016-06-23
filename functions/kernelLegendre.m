function [S]=kernelLegendre(Nmax,L)
%{
 Recive:
 Nmax: Orden maximo de polinomio que se genera, 
 L:    Longitud de los vectores

 Entrega: S Matriz (Nmax+1 x L)
%}

%format long

% Coefficients
Aii=zeros(Nmax+1,Nmax+2);
Z=zeros(Nmax-2,Nmax-1);
Aii(1,2)=1;
Aii(2,3)=1;
for k=2:Nmax
    f1=(2*k-1)/k;
    f2=-(k-1)/k;
    for i=2:k+2
        Aii(k+1,i)=f1*Aii(k,i-1)+f2*Aii(k-1,i);
    end
end
Aii=Aii(1:end,2:end);

% Evaluation

x=linspace(-1,1,L);
% x=horzcat(x,0);
S=zeros((Nmax+1),length(x));

for k=1:(Nmax+1)
    y=polyval(Aii(k,k:-1:1),x);
    Ynorm=norm(y);
    y=y/Ynorm;
    S(k,:)=y;
end

end

        