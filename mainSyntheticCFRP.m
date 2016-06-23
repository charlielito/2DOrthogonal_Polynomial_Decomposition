 close all; clear all; clc;

% Add to path 'functions' directory from whereever you run this 
% (obviosly that folder must be with this file)
currdir = mfilename('fullpath');
sep = strfind(currdir,'\');
path(path,[currdir(1:sep(end)),'functions']);


disp('-------- Loading data... ------')
load syntheticCFRP %RT and t and RT_coeffall and n
disp('-------- Data loaded! ------')
disp('-------- Calculating things... ------')

mc = [250 50  50  250 250 50 50  250 250]; %1 pixel approx 1mm
nc = [150 200 100 250 50  50 250 100 200]; 
Tmat = RT;
P = size(RT,3);
N = length(mc);
fs = 1/0.018;
logt = log((1:P)/fs);


ls = 15*1.5; %long side of window

%For POD windows
mi = mc - round(ls/2);
mf = mc + round(ls/2);
ni = nc - round(ls/2);
nf = nc + round(ls/2);
K0mn = zeros(1,P);

ls = round(ls);

figure(1),imagesc(RT(:,:,400))

ms = 234;
ns = 170;
Ts = zeros(P,1);
win = 5;
hold on
for j = 1:win
    for i = 1:win
        Ts = Ts + squeeze(RT(ms+j-1,ns+i-1,:));
        plot(ns+j-1,ms+i-1,'.k','MarkerSize',3)
    end
end
hold off
axis off

Ts = Ts/(win*win);

Ts0 = Ts/Ts(1);

times = zeros(1,N);
times2 = times;
times3 = times;
depths = 0.2:0.2:1.8;


for i = 1:N
    figure(2)
    subplot(2,5,i)
    Tdef = squeeze(Tmat(mc(i), nc(i),:))./squeeze(Tmat(mc(i), nc(i),1)) - Ts0;
    
    temp = find(Tdef == max(Tdef));
    if temp/fs < exp(logt(end-10)) %there is no right contrast peak
    times2(i) = temp; end
    plot(logt,Tdef,'k',logt(temp), Tdef(temp),'xr')
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlabel('ln(t)')
    ylabel('C^n(t)')
    xlim([min(logt) max(logt) ])
    %ylim([max(Tdef)*0.5 max(Tdef)*1.1 ])

    
    %Calculate 1st derivatives
    n = 1;
    logt = log((1:P)/fs);
    logt = logt(n:end);
    Tdef = log(squeeze(Tmat(mc(i), nc(i),n:end)));
    Ts = log(Ts0(n:end));
    p = polyfit(logt',Tdef,7);
    Tdefp = polyval(p,logt);
    p = diffPoly(p);
    DTdef = polyval(p,logt);
    ps = polyfit(logt',Ts,7);
    ps = diffPoly(ps);
    DTs = polyval(ps,logt);

    %Calculates 2nd derivatives
    figure(3)
    subplot(2,5,i)
    p = diffPoly(p);
    DDTdef = polyval(p,logt);
    ps = diffPoly(ps);
    DDTs = polyval(ps,logt);
    [Max, Idx] = findpeaks(DDTdef);
   
    if ~isempty(Idx) times(i) = Idx(1);
    else
    Idx = 1;Max = DDTdef(Idx); end
    plot(logt,DDTdef,'k',logt(Idx), Max, 'xr')
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlabel('ln(t)')
    ylabel('d^2ln(T)/d[ln(t)] ^2')
    xlim([min(logt) max(logt) ])
    ylim([-0.05 0.05])

    % Calculate OPD
    logt = log((1:P)/fs);
    figure(4)
    subplot(2,5,i)
    intM = mi(i):mf(i);
    intN = ni(i):nf(i);

    
    for pp = 1:P
        Matrix_mn = squeeze(RT(intM,intN,pp));
        Coeffmn = OPD((Matrix_mn),[ls ls]);
        K0mn(pp) = (mean(abs(Coeffmn(3,:))) + mean(abs(Coeffmn(:,3))))/2;

    end
    [fpol, fgauss, tMin, tMax ] = remove_tendency(K0mn,10);
    plot(logt,log(K0mn)-min(log(K0mn)),'k',logt,fgauss,'-.k',logt(tMax),fgauss(tMax),'xr')
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlabel('ln(t)')
    ylabel('ln(Coeff)')
    xlim([min(logt) max(logt) ])
    times3(i) = tMax;

end

%%
step = 0.018;
limit = 5;
depths = 0.2:0.2:1.8;
t1 = times(1:limit)*step;
figure, plot(t1, depths(1:limit),'xr', 'linewidth',1.5,'MarkerSize',10)
set(get(gca,'xlabel'),'string','Specific characteristic time (s)','fontsize',12)
set(get(gca,'ylabel'),'string','Depth (mm)','fontsize',12)
set(get(gca,'title'),'string','Depth estimation','fontsize',11, 'FontWeight','bold')

degree = 2;
[p,s] = polyfit(t1', depths(1:limit)',degree);
t2 = linspace(3.5,11.5);
z = polyval(p,t2);
hold on
h1 = plot(t2,z,'r');
R = zeros(1,3);
R(1) =  (1 - s.normr^2/norm(depths(1:limit)-mean(depths(1:limit)))^2)*100;
Zestimated = zeros(3,length(depths));
Zestimated(1,:) = polyval(p,times*step);
p_tsr = p;


t1 = times2*step;
plot(t1, depths,'ob', 'linewidth',1.5,'MarkerSize',7)
[p,s] = polyfit(t1(~isnan(t1))', depths(~isnan(t1))',1);
t2 = linspace(0,10.5);
z = polyval(p,t2);
hold on
h2 = plot(t2,z,'.-.b','linewidth',1.5,'markersize',3);
R(2) =  (1 - s.normr^2/norm(depths-mean(depths))^2)*100;
Zestimated(2,:) = polyval(p,t1);
p_ntc = p;


linear = 0; %linear fit or not linear (sqrt(x))
t1 = times3*step;
plot(t1, depths,'+k', 'linewidth',1.5,'MarkerSize',7)
t2 = linspace(0.4,10);
if (linear)
    [p,s] = polyfit(t1', depths',1);
    z = polyval(p,t2);
    R(3) =  (1 - s.normr^2/norm(depths-mean(depths))^2)*100;
    Zestimated(3,:) = polyval(p,t1);
else
    Eqn = 'a*sqrt(x) + b';
    [f1, gof] = fit(t1',depths',Eqn,'Startpoint',[1 1]);
    z = f1.a*sqrt(t2)+f1.b;
    R(3) = gof.rsquare*100;
    Zestimated(3,:) = f1.a*sqrt(t1)+f1.b;
end
hold on
h3 = plot(t2,z,':k','linewidth',2);


text(t1(1),depths(end)-0.4,['R^2 = ', num2str(R(3)), ' %'],'Color','black')
text(t1(1),depths(end)-0.2,['R^2 = ', num2str(R(2)), ' %'],'Color','blue')
text(t1(1),depths(end),['R^2 = ', num2str(R(1)), ' %'],'Color','red')
legend([h1 h2 h3],'TSR', 'NTC', 'OPD','Location','southeast')


Zestimated;
Zexact = [depths;depths;depths];
Zerror = (100*abs(Zestimated-Zexact)./Zexact);
clc
disp('*********** Predicted values (mm) *******')
disp('Value     TSR         NTC         OPD')
disp(num2str([depths' Zestimated']))
disp('*********** Relative errors (%) *******')
disp('  TSR           NTC           OPD')
disp(num2str(Zerror'))