% Transforms a 3D matrix MxNxP into a 2D matrix (M*N)xP matrix
% The columns of the new matrix contain each "frame" of the 3D matrix
% a "standardisation" process is applied to each column vector: 
 %      each pixel A(n,m)  = (A(n,m) - u)/standardev

function newM = compress3DTo2DForPCA(data)
    [M, N, P] = size(data);
    newM = zeros(M*N,P);
    for i = 1:P
        aux = data(:,:,i);
        aux =  aux'; % to expand from row by row
        colvect = aux(:); %apply standardisation
        newM(:,i) = (colvect-mean(colvect))/std(colvect);
    end
    
end



