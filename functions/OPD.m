
function C = OPD(mat,P, Phi_m, Phi_n)

       
[Ml,Nl] = size(mat);

Mp = P(1);
Np = P(2);

if (nargin == 2)
    % Kernel for basis functions
    [Phi_m]=kernelLegendre(Ml-1,Ml);
    [Phi_n]=kernelLegendre(Nl-1,Nl);

    [Phi_m, Crr]=qr(Phi_m');
    [Phi_n, Crr]=qr(Phi_n');

    Phi_m=Phi_m(1:end,1:(Mp))';
    Phi_n=Phi_n(1:end,1:(Np))';
end

% Matrix with coefficients             
C = Phi_m*mat*Phi_n';  