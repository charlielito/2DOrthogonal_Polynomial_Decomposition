
function [fpol, fgauss, tMin, tMax ] = remove_tendency(mat,degree)

N = length(mat);
t = log(1:N);
logK = log(mat);

%find the minimun to shift all to cero DC lvl.
globalmin = min(logK);
y = logK;%-globalmin;

%polynomial fit to find local minima
[p,s,u] = polyfit(t,y,degree); %10th order works fine!
fpol = polyval(p,t,[],u);

nopoint = 0;
DataInv = 1.01*max(fpol) - fpol;
[~,MinIdx] = findpeaks(DataInv);

if (~isempty(MinIdx)) 
maxv = MinIdx(1);
numax = 1;
if maxv > N/5 maxv = numax;end
else
    nopoint = 1;
    maxv = 1;
end
%maxv


t = t(1:maxv);
logK2 = logK(1:maxv);
y = logK2-globalmin;

% and one last point of tail
% if (~nopoint)
% npoints = 100;
% y = [y logK(end-npoints:end)-globalmin];
% t = [t log(N-npoints:N)];    
% else
% y = [y mean(logK(end-150:end))-globalmin];
% t = [t log(N)];
% end

npoints = 20;
y = [y logK(end-npoints:end)-globalmin];
t = [t log(N-npoints:N)];  

%Exponential fit
ExpEqn = 'a*exp(-b*x)+c'; 
f1 = fit(t',y',ExpEqn, 'Lower',[0,0,-Inf],'Upper',[Inf,Inf,Inf],'Startpoint',[1 1 1]);


%Now with all data
n = 1:N;
t = log(n);
y = logK-globalmin;

%Subtract Noise
NoiseModel = f1.a*exp(-f1.b*t)+f1.c;
newy = y-NoiseModel;


gaussEqn = 'a*exp(-((x-b)/c)^2)+d'; %therefore b is mean, c = 4sigma, a peak value
indexMax = find( newy == max(newy));
logTindexMax = t(indexMax); %in case there are two
% f2 = fit(t',newy',gaussEqn, 'Lower',[0,0,-Inf,-Inf],'Upper',[Inf,max(t),max(t),Inf],...
%     'Startpoint',[1 t(indexMax) 1 1]);

f2 = fit(t(maxv:end)',newy(maxv:end)',gaussEqn, 'Lower',[0,0,-Inf,-Inf],'Upper',[Inf,max(t),max(t),Inf],...
    'Startpoint',[1 logTindexMax(1) 1 1]);


tMin = (maxv);
tMax = round(exp(f2.b));
fgauss = f2.a*exp(-((t-f2.b)/f2.c).^2)+f2.d;
%figure,plot(t,y,t,newy,t,fgauss)
end
