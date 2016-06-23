function pol = diffPoly(p)

    N =  length(p);
    pol = zeros(1,N-1);
    for i = 1:N-1
        pol(i) = (N-i)*p(i);
    end
end