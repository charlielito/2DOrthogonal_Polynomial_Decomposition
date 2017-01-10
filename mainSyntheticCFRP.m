 close all; clear all; clc;

% Add to path 'functions' directory from whereever you run this 
% (obviosly that folder must be with this file)
currdir = mfilename('fullpath');
sep = strfind(currdir,'\');
path(path,[currdir(1:sep(end)),'functions']);


disp('-------- Loading data... ------')
load syntheticCFRP %RT and t
disp('-------- Data loaded! ------')
disp('-------- Calculating things... ------')

mc = [250 50  50  250 250 50 50  250 250]; %1 pixel approx 1mm
nc = [150 200 100 250 50  50 250 100 200]; 
P = size(RT,3);
N = length(mc);
fs = 1/0.018;
logt = log((1:P)/fs);


ls = 15*2; %long side of window

%For POD windows
mi = mc - round(ls/2);
mf = mc + round(ls/2);
ni = nc - round(ls/2);
nf = nc + round(ls/2);
K0mn = zeros(1,P);

ls = round(ls);

figure(1),imagesc(RT(:,:,150))

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
axis square

Ts = Ts/(win*win);

Ts0 = Ts/Ts(1);

times = zeros(1,N);
times2 = times;
times3 = times;
times4 = times;
times5 = times;
times3_2 = times;
depths = 0.2:0.2:1.8;

%PPT parameters
[rows, cols, frames] = size(RT);
Fs = 55;
nppt = 1;
Tsppt = Ts;

for i = 1:N
    % ****** NTC ***********
    figure(2)
    subplot(2,5,i)
    Tdef = squeeze(RT(mc(i), nc(i),:))./squeeze(RT(mc(i), nc(i),1)) - Ts0;
    
    temp = find(Tdef == max(Tdef));
    if temp/fs < exp(logt(end-10)) %there is no right contrast peak
    times2(i) = temp; end
    plot(logt,Tdef,'k',logt(temp), Tdef(temp),'xr')
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlabel('ln(t)')
    ylabel('C^n(t)')
    xlim([min(logt) max(logt) ])

    % ****** TSR ***********
    %Calculate 1st derivatives
    n = 1;
    logt = log((1:P)/fs);
    logt = logt(n:end);
    Tdef = log(squeeze(RT(mc(i), nc(i),n:end)));
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
    [Max, Idx] = findpeaks(DDTdef,'SORTSTR','ascend');
   
    if ~isempty(Idx) times(i) = Idx(1);
    else
    Idx = 1;Max = DDTdef(Idx); end
    plot(logt,DDTdef,'k',logt(Idx), Max, 'xr'), hold all
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlabel('ln(t)')
    ylabel('d^2ln(T)/d[ln(t)] ^2')
    xlim([min(logt) max(logt) ])
    ylim([-0.05 0.09])

    % --------------- Calculate OPD and PCT
    intM = mi(i):mf(i);
    intN = ni(i):nf(i);
    logt = log((1:P)/fs);
    
    % *****PCT **********
    figure(4)
    subplot(2,5,i)
    Matrix_mn = RT(intM,intN,:);
    PCA_Mat = compress3DTo2DForPCA(Matrix_mn);
    [U,S,V] = svd(PCA_Mat);
    x = V(:,2)'; %time profile from 2nd orthogonal basis
    [Max, Idx] = findpeaks(x,'SORTSTR','ascend','minpeakdistance',100);
    DataInv = 1.01*max(x) - x;
    [Minima,MinIdx] = findpeaks(DataInv,'SORTSTR','ascend','minpeakdistance',100);
    plot(logt,x,'k')
    xlabel('ln(t)')
    ylabel('2nd Orthogonal Vector')
    hold on
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlim([min(logt) max(logt) ])
    if isempty(Max) || Idx(end) < 15 || Idx(end) > 300% no maxima or fake maxima
        times4(i) = MinIdx(end);
        plot(logt(MinIdx(end)), x(MinIdx(end)),'rx')
    else
        times4(i) = Idx(end);
        plot(logt(Idx(end)),Max(end),'rx')
    end
    hold off

    % **** OPD ************
    figure(5)
    subplot(2,5,i)
    for pp = 1:P
        Matrix_mn = squeeze(RT(intM,intN,pp));
        Coeffmn = OPD((Matrix_mn),[ls ls]);
        K0mn(pp) = (mean(abs(Coeffmn(3,:))) + mean(abs(Coeffmn(:,3))))/2;

    end
    [fpol, fgauss, tMin, tMax ] = remove_tendency(K0mn,12);
    plot(logt,log(K0mn)-min(log(K0mn)),'k',logt,fgauss,'-.k',logt(tMax),fgauss(tMax),'xr')
    %plot(logt,fgauss)
    hold all
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlabel('ln(t)')
    ylabel('ln(Coeff)')
    xlim([min(logt) max(logt) ])
    times3(i) = tMax;
    [~, Index] = findpeaks(fpol,'SORTSTR','ascend');
    if ~isempty(Index)    times3_2(i) = Index(end); end

    
% ******** PPT ********************
    figure(6)
    Temp = squeeze(RT(mc(i), nc(i),1:nppt:frames));
    TdefFreq = fft(Temp);
    TsFreq = fft(Tsppt(1:nppt:frames));
    Frames = size(TdefFreq,1);
    TdefFreq = TdefFreq(1:floor(Frames/2)+1);%only positive frequencies
    TsFreq = TsFreq(1:floor(Frames/2)+1);%only positive frequencies
    f = (0:Frames/2)*Fs/Frames/nppt;
    Pdef = angle(TdefFreq); %phase profiles
    Psound = angle(TsFreq);
    PCon = Psound-Pdef; %phase contrast
    factor = 10; %resampling factor to have better blind frequency estimation
    PCon2 = resample(PCon, factor,1);
    f2 = (0:length(PCon2)-1)*Fs/Frames/nppt/factor; %frequency resampled
    times5(i) = f2(zerocrossPPT(PCon2(2:end))); %calculate zero cross (blind f_b)
    plot(f2,PCon2)
    xlabel('Frequency (Hz)')
    ylabel('Phase Contrast (rad)')
    hold all
    title('Phase contrasts for different dephts ')
    xlim([min(f) max(f) ])
end
figure(6),legend('0.2 mm', '0.4 mm', '0.6 mm','0.8 mm','1.0 mm','1.0 mm','1.2 mm','1.4 mm','1.6 mm','1.8 mm')
times3_2(8) = exp(0.775)/0.018;

%% ****** TSR ***********
step = 0.018;
limit = 4;
depths = 0.2:0.2:1.8;
t1 = times(1:limit)*step;
figure, plot(t1, depths(1:limit),'xr', 'linewidth',1.5,'MarkerSize',10)
set(get(gca,'xlabel'),'string','Specific characteristic time t (s)','fontsize',12)
set(get(gca,'ylabel'),'string','Depth z (mm)','fontsize',12)
set(get(gca,'title'),'string','Depth estimation','fontsize',11, 'FontWeight','bold')

degree = 2;
[p,s] = polyfit(t1', depths(1:limit)',degree);
%t2 = linspace(3.5,11.5);
t2 = linspace(1,9);
z = polyval(p,t2);
hold on
h1 = plot(t2,z,'r','linewidth',1.5,'markersize',3);
R = zeros(1,5);
R(1) =  (1 - s.normr^2/norm(depths(1:limit)-mean(depths(1:limit)))^2)*100;
Zestimated = zeros(5,length(depths));
Zestimated(1,:) = polyval(p,times*step);
p_tsr = p;


% ****** NTC ***********
limit = 7;
depths = depths(1:limit);
t1 = times2(1:limit)*step;
plot(t1, depths,'ob', 'linewidth',1.5,'MarkerSize',7)
[p,s] = polyfit(t1(~isnan(t1))', depths(~isnan(t1))',1);
t2 = linspace(0,10.5);
z = polyval(p,t2);
hold on
h2 = plot(t2,z,'.-.b','linewidth',1.5,'markersize',3);
R(2) =  (1 - s.normr^2/norm(depths-mean(depths))^2)*100;
Zestimated(2,:) = polyval(p,times2*step);
p_ntc = p;
depths = 0.2:0.2:1.8;

% ****** OPD ***********
%limit = 8;
depths = depths(1:limit);
linear = 1; %linear fit or not linear (sqrt(x))
t1 = times3(1:limit)*step;
plot(t1, depths,'+k', 'linewidth',1.5,'MarkerSize',7)
t2 = linspace(0.0,11);
if (linear)
    [p,s] = polyfit(t1', depths',1);
    z = polyval(p,t2);
    R(3) =  (1 - s.normr^2/norm(depths-mean(depths))^2)*100;
    Zestimated(3,:) = polyval(p,times3*step);
    p_opd = p;
else
    Eqn = 'a*sqrt(x) + b';
    [f1, gof] = fit(t1',depths',Eqn,'Startpoint',[1 1]);
    z = f1.a*sqrt(t2)+f1.b;
    R(3) = gof.rsquare*100;
    Zestimated(3,:) = f1.a*sqrt(times3*step)+f1.b;
    p_opd(1) = f1.a;
    p_opd(2) = f1.b;
end
hold on
h3 = plot(t2,z,':k','linewidth',2);
depths = 0.2:0.2:1.8;


% ****** PCT ***********
%limit = 8;
depths = depths(1:limit);
linear = 1; %linear fit or not linear (sqrt(x))
t1 = times4(1:limit)*step;
plot(t1, depths,'*g', 'linewidth',1.5,'MarkerSize',7)
t2 = linspace(0,11);
if (linear)
    [p,s] = polyfit(t1', depths',1);
    z = polyval(p,t2);
    R(4) =  (1 - s.normr^2/norm(depths-mean(depths))^2)*100;
    Zestimated(4,:) = polyval(p,times4*step);
    p_pca = p;
else
    Eqn = 'a*sqrt(x) + b';
    [f2, gof] = fit(t1',depths',Eqn,'Startpoint',[1 1]);
    z = f2.a*sqrt(t2)+f2.b;
    R(4) = gof.rsquare*100;
    Zestimated(4,:) = f2.a*sqrt(times4*step)+f2.b;
    p_pca(1) = f2.a;
    p_pca(2) = f2.b;
end
hold on
h4 = plot(t2,z,'--g','linewidth',1.5,'markersize',3);
depths = 0.2:0.2:1.8;


% **** PPT *************
range = 1:7;
t = times5(range).^(-1);
plot(t, depths(range),'^m', 'linewidth',1.5,'MarkerSize',7)
[f1,gof] = fit(t',depths(range)','a*sqrt(x) + b');
t2 = linspace(0,6);
z = f1.a*sqrt(t2) + f1.b;
h5 = plot(t2,z,'m');
R(5) =  gof.rsquare*100;
p_ppt = [f1.a f1.b];
Zestimated(5,1:length(times5)) = p_ppt(1)*sqrt(times5.^-1) + p_ppt(2);

%round to 1 and 2 decimals
R = (round(10*R)/10);
p_pca = (round(1000*p_pca)/1000);
p_tsr = (round(1000*p_tsr)/1000);
p_ntc = (round(1000*p_ntc)/1000);
p_opd = (round(1000*p_opd)/1000);
p_ppt = (round(1000*p_ppt)/1000);


% ----Place equations and errors R^2
depths = 0.2:0.2:1.8;
signn = containers.Map({1, -1},{'+','-'});

t1(1)=t1(1)-0.2;
%TSR
text(t1(1),depths(end)+0.1,['R^2_{TSR} = ', num2str(R(1)), ' %'],'Color','red','FontSize',8)
text(t1(1),depths(end),['z_{TSR} = ', num2str(p_tsr(1)), 't^2 ', ...
    signn(sign(p_tsr(2))),' ', num2str(abs(p_tsr(2))), 't ', ...
    signn(sign(p_tsr(3))),' ',  num2str(abs(p_tsr(3))) ],'Color','red','FontSize',8)
%NTC
text(t1(1),depths(end)-0.15,['R^2_{NTC} = ', num2str(R(2)), ' %'],'Color','blue','FontSize',8)
text(t1(1),depths(end)-0.25,['z_{NTC} = ', num2str(p_ntc(1)), 't ', ...
    signn(sign(p_ntc(2))),' ',  num2str(abs(p_ntc(2))) ],'Color','blue','FontSize',8)
%OPD
text(t1(1),depths(end)-0.4,['R^2_{OPD} = ', num2str(R(3)), ' %'],'Color','black','FontSize',8)
text(t1(1),depths(end)-0.5,['z_{OPD} = ', num2str(p_opd(1)), 't ', ...
    signn(sign(p_opd(2))),' ',  num2str(abs(p_opd(2))) ],'Color','black','FontSize',8)
%PCT
text(t1(1),depths(end)-0.65,['R^2_{PCT} = ', num2str(R(4)), ' %'],'Color','green','FontSize',8)
text(t1(1),depths(end)-0.75,['z_{PCT} = ', num2str(p_pca(1)), 't ', ...
    signn(sign(p_pca(2))),' ',  num2str(abs(p_pca(2))) ],'Color','green','FontSize',8)
% PPT
text(t1(1),depths(end)-0.9,['R^2_{PPT} = ', num2str(R(5)), ' %'],'Color','m','FontSize',8)
text(t1(1),depths(end)-1.0,['z_{PPT} = ', num2str(p_ppt(1)), 't^{1/2} '],'Color','m','FontSize',8)

legend([h1 h2 h3 h4 h5],'TSR', 'NTC', 'OPD','PCT', 'PPT','Location','southeast')
axis( [0 6 0.1 1.99])

% Calculate Errors and print them
Zestimated;
Zexact = [depths;depths;depths;depths;depths];
Zerror = (100*(round(100*Zestimated)/100-Zexact)./Zexact);
clc
disp('*********** Predicted values (mm) *******')
disp('Value        TSR          NTC          OPD          PCT          PPT')
disp(num2str([depths' Zestimated'],'%0.2f        '))
disp('*********** Relative errors (%) *********')
disp('             TSR          NTC          OPD          PCT          PPT')
disp(num2str([nan*ones(1,size(Zerror,2))' Zerror'],'%0.1f       '))