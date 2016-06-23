function new_sequence = ntc(T,t,winsize,x,y,mc,nc)

% T : 3D array which corresponds to a thermograms secuence

for i = 1:winsize;
    for j = 1:winsize
        di((i-1)*winsize+j) = x+j-1;
        dj((i-1)*winsize+j) = y+i-1;
    end
end
i = winsize*winsize;

begin = 2;
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

tm  = 1;
Ts = pixelsequence/pixelsequence(tm);
M_mc = size(mc,1);
N_mc = size(mc,2);
new_sequence = zeros(M_mc,N_mc,n_f);
for mm = 1:M_mc
    for nn = 1:N_mc
        y = log(T(mc(mm,nn),nc(mm,nn),:));
        ypol = polynomialfit(logt(:),y(:),n);
        Tdef = polyval(ypol,logt);
        Tdef = exp(Tdef(begin:end)); %back to Temperature
        new_sequence(mm,nn,:) = Tdef/Tdef(tm) - Ts;
    end
end
        


end