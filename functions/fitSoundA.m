function pixelsequence = fitSoundA(T,t,winsize,x,y)

% T : 3D array which corresponds to a thermograms secuence

for i = 1:winsize;
    for j = 1:winsize
        di((i-1)*winsize+j) = x+j-1;
        dj((i-1)*winsize+j) = y+i-1;
    end
end
i = winsize*winsize;

begin = 1;
[N,l,n_f] = size(T(:,:,begin:end));
pixelsequence = zeros(1,n_f);
n = 9;
logt = log(t);
for j = 1 : i-1
    y = log(T(dj(j),di(j),:));
    ypol = polynomialfit(logt(:),y(:),n);
    yfit = polyval(ypol,logt);
    pixelsequence = pixelsequence + exp(yfit(begin:end));
end
pixelsequence = pixelsequence/j;

end