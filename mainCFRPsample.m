close all; clear all; clc;
% Add to path 'functions' directory from whereever you run this 
% (obviosly that folder must be with this file)
currdir = mfilename('fullpath');
sep = strfind(currdir,'\');
path(path,[currdir(1:sep(end)),'functions']);



disp('-------- Loading data... ------')
load CFRPsample RTfront tfront RTrear trear
disp('-------- Data loaded! ------')

[m,n,p] = size(RTfront);
[m2,n2,p2] = size(RTrear);
p = min([p p2]);

% Front Sample defects Centers
mc = [73   152   228   304   380];
nc = [85   164   238   314   394];
%to fit better centers due to deformities of material
offsetm  = [0 0 0 3 0; 2 -3 0 0 2; 2 0 0 0 2; 4 -2 2 0 4; 8 0 0 0 8];
offsetn  = [0 0 0 2 0; 2 0 0 2 2; 0 0 0 0 0; 0 0 0 0 0; -2 0 2 0 0];

dm = mc(2:end) - mc(1:(end-1));
dm = mean(dm);
dn = nc(2:end) - nc(1:(end-1));
dn = mean(dn);
d = mean([dm dn]);

ls = [3 15 3 15 3;
      5 10 5 10 5;
      7  7 7  7 7;
      10 5 10 5 10;
      15 3 15 3 15]; %sizes: 3 smallest one
                       %    15 biggest one                     
kk = 1.5;
mls = kk*round(ls*d/50);
pmn = round(mls/2); %polynomial order
%pmn = mls; 
mc = [mc' mc' mc' mc' mc'];
nc = [nc;nc;nc;nc;nc];
mc = mc + offsetm;
nc = nc + offsetn;
mi = mc - round(0.5*mls);
mf = mc + round(0.5*mls);
ni = nc - round(0.5*mls);
nf = nc + round(0.5*mls);


% Back Sample defects Centers
mc2 = [46   130   207   286   370];
nc2 = [70   150   238   321   394];
%to fit better centers due to deformities of material
offsetm  = [0 0 0 0 0; -2 0 0 0 2; 2 0 0 0 5; 0 0 0 0 5; -2 0 0 0 2];
offsetn  = [0 0 0 0 0; 2 0 0 0 -1; 1 0 0 0 -1; 2 0 0 0 0; 0 0 0 0 0];

dm = mc2(2:end) - mc2(1:(end-1));
dm = mean(dm);
dn = nc2(2:end) - nc2(1:(end-1));
dn = mean(dn);
d = mean([dm dn]);

mls = kk*round(ls*d/50);
pmn2 = round(mls/2); %polynomial order
mc2 = [mc2' mc2' mc2' mc2' mc2'];
nc2 = [nc2;nc2;nc2;nc2;nc2];
mc2 = mc2 + offsetm;
nc2 = nc2 + offsetn;
mi2 = mc2 - round(0.5*mls);
mf2 = mc2 + round(0.5*mls);
ni2 = nc2 - round(0.5*mls);
nf2 = nc2 + round(0.5*mls);


 
K0mn = zeros(5,5,p);
K0mn2 = zeros(5,5,p);

disp('-------- Calculating OPD for all defects... ------')

for pp = 1:p

for mm = 1:5
    for nn = 1:5
        intM = mi(mm,nn):mf(mm,nn);
        intN = ni(mm,nn):nf(mm,nn);
        intM2 = mi2(mm,nn):mf2(mm,nn);
        intN2 = ni2(mm,nn):nf2(mm,nn);

        Matrix_mn = squeeze(RTfront(intM,intN,1+(pp-1)));
        Matrix_mn2 = squeeze(RTrear(intM2,intN2,1+(pp-1)));

        Coeffmn = OPD((Matrix_mn),[pmn(mm,nn) pmn(mm,nn)]);
        Coeffmn2 = OPD((Matrix_mn2),[pmn2(mm,nn) pmn2(mm,nn)]);
        K0mn(mm,nn,pp) = (mean(abs(Coeffmn(3,:))) + mean(abs(Coeffmn(:,3))))/2;
        K0mn2(mm,nn,pp) = (mean(abs(Coeffmn2(3,:))) + mean(abs(Coeffmn2(:,3))))/2;
    end
end

end


disp('-------- OPD calculated! ------')
%%
disp('-------- Calculating Normalized Contrast... ------')
NCfront = ntc(RTfront,tfront,5,100,380,mc,nc);
disp('-------- Normalized Contrast for front part terminated! ------')
NC = ntc(RTrear,trear,5,90,360,mc2,nc2); %for rear sample
disp('-------- Normalized Contrast for REAR part terminated! ------')

%%
fs = 1/55; %sample period 

Sm = [1 5 1 5 1;   % Position according to sizes
      2 4 2 4 2;
      3 3 3 3 3;
      4 2 4 2 4;
      5 1 5 1 5];
  
Hn = [3 4 2 5 1];  % Depths for front
Hn2 = [5 1 4 2 3];  % Depths for rear

% mm indicates size, 1 -> smallest, 5-> biggest
% nn indicates depth, 1-> shallowest, 5-> deepest
% Acces the is -> K( Sm(mm,Hn(nn)) ,Hn(nn) )

mm = 5; %biigest flaw
CoeffsFront = [];
CoeffsBack = [];
NCbigFront = [];
NCbigBack = [];
RTbigFront = [];
RTbigBack = [];
for nn = 1:5
    %OPD
    CoeffsFront = [CoeffsFront; squeeze(K0mn( Sm(mm,Hn(nn)), Hn(nn), :))'];
    CoeffsBack = [CoeffsBack; squeeze(K0mn2( Sm(mm,Hn2(nn)), Hn2(nn), :))'];
    
    %NTC
    NCbigFront = [NCbigFront; squeeze(NCfront( Sm(mm,Hn(nn)), Hn(nn), :))'];
    NCbigBack = [NCbigBack; squeeze(NC( Sm(mm,Hn2(nn)), Hn2(nn), :))'];

    %TSR
    RTbigFront = [RTbigFront; squeeze(RTfront( mc(Sm(mm,Hn(nn)), Hn(nn)), nc(Sm(mm,Hn(nn)), Hn(nn)), :))'];
    RTbigBack = [RTbigBack; squeeze(RTrear( mc2(Sm(mm,Hn2(nn)), Hn2(nn)), nc2(Sm(mm,Hn2(nn)), Hn2(nn)), :))'];
end

N = min([length(CoeffsFront), length(CoeffsBack)]);
CoefTotal = [CoeffsFront(:,1:N); CoeffsBack(:,1:N)];

Nn = min([length(NCbigFront), length(NCbigBack)]);
NCbig = [NCbigFront(:,1:Nn); NCbigBack(:,1:Nn)];


Nn = min([length(RTbigFront), length(RTbigBack)]);
RTbig = [RTbigFront(:,1:Nn); RTbigBack(:,1:Nn)];

depths = [0.2:0.2:1 1:0.2:1.8];

N = length(CoefTotal(1,:));
t5 = log((1:N)*fs);
%t = log10(1:N);
Ndefect = 10;
tGauss = zeros(1,Ndefect);
tPol = zeros(1,Ndefect);

f = size(RTbig,2);
t = log((1:f)*fs);
t2 = 1:f;
t3 = log(1:f);
times = zeros(1,10);
times2 = zeros(1,10);
for i = 1:Ndefect
    figure(1)
    subplot(2,5,i)
    C = CoefTotal(i,:);
    [fpol, fgauss, tMin, tMax ] = remove_tendency(C,10); 
    [~,MinIdx] = findpeaks(fpol);
    if (length(MinIdx) >= 3) tminp = MinIdx(2);
    elseif (isempty(MinIdx)) tminp = nan;
    else tminp = min(MinIdx);end
    tPol(i) = tminp*fs;
    tGauss(i) = tMax*fs;
    plot(t5,log10(C),t5,fpol,'r',t5,fgauss,'g')
    plot(t5,log(C)-min(log(C)),'k',t5,fgauss,'-.k',t5(tMax),fgauss(tMax),'xr')
    xlabel('ln(t)')
    ylabel('ln(Coeff)')
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlim([min(t5) max(t5) ])
    
    
    figure(2)
    subplot(2,5,i)
    
    ncc = NCbig(i,:);
    f1 = polyfit(t(1:end-1)',ncc(1:end)',7);
    pNC = polyval(f1,t);
    [Max, Indx] = findpeaks(double(pNC),'NPEAKS',10);

    tmax = Indx( Max == max(Max(Max<Inf) ) );
    ear = max(ncc)*0.1+min(ncc);
    timeSI = Indx( Max == max(Max(Max<Inf) ) );
    times2(i) = fs*timeSI;

    plot(t,pNC,'k',t(timeSI),pNC(timeSI),'xr')
    xlabel('ln(t)')
    ylabel('C^n(t)')
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlim([min(t) max(t) ])
   

    Tdef = log(RTbig(i,2:end));
    t3 = log((2:f)*fs);
    t3 = t3(1:end);
    p = polyfit(t3',Tdef',8);
    p = diffPoly(p);
    DTdef = polyval(p,t3);

    
    figure(3)
    subplot(2,5,i)
    p = diffPoly(p);
    DDTdef = polyval(p,t3);
    [Max, Indx] = findpeaks(double(DDTdef));
    if length(Indx) == 1
        times(i) = Indx(1)*fs;
    else
    times(i) = Indx(2)*fs;
    end
    plot(t3,DDTdef,'k',t3(times(i)/fs),DDTdef(times(i)/fs),'xr');
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlabel('ln(t)')
    ylabel('d^2ln(T)/d[ln(t)] ^2')
    xlim([min(log((1:f)*fs)) max(t3) ])
end
%%
depths = 0.2:0.2:1.8;
limit = 4;
t1 = times(1:limit);
figure(4), plot(t1, depths(1:limit),'xr', 'linewidth',1.5,'MarkerSize',10), hold on

degree = 2;
tf = 6;
[p,s] = polyfit(t1', depths(1:limit)',degree);
t2 = linspace(t1(1),tf);
z = polyval(p,t2);
h1 = plot(t2,z,'r');
R = zeros(1,3);
R(1) =  (1 - s.normr^2/norm(depths(1:limit)-mean(depths(1:limit)))^2)*100;
Zestimated = zeros(3,length(depths));
Zestimated(1,1:length(t1)) = polyval(p,t1);
p_tsr = p;


tf = 2.5;
Ndefect = 7;
degree = 1;
depths = 0.2:0.2:0.2*(Ndefect-1);
times22 = [times2(1:4), mean(times2(6:6)),times2(7:Ndefect)];
plot(times22, depths,'ob', 'linewidth',1.5,'MarkerSize',7)
[p,s] = polyfit(times22', depths',degree);
t2 = linspace(0,tf);
z = polyval(p,t2);
h2 = plot(t2,z,'.-.b','linewidth',1.5,'markersize',3);
R(2) =  (1 - s.normr^2/norm(depths-mean(depths))^2)*100;
Zestimated(2,1:length(times22)) = polyval(p,times22);
p_ntc = p;


depths = 0.2:0.2:1.2;
tGauss2 = [tGauss(1:4), tGauss(6:7)];
plot(tGauss2, depths,'+k', 'linewidth',1.5,'MarkerSize',8)
t1 = tGauss2;
t2 = linspace(0,2.5);
linear = 1;
if (linear)
    [p,s] = polyfit(t1', depths',1);
    z = polyval(p,t2);
    R(3) =  (1 - s.normr^2/norm(depths-mean(depths))^2)*100;
    Zestimated(3,1:length(t1)) = polyval(p,t1);
else
    Eqn = 'a*sqrt(x) + b';
    [f1, gof] = fit(t1',depths',Eqn);
    z = f1.a*sqrt(t2)+f1.b;
    R(3) = gof.rsquare*100;
    Zestimated(3,1:length(t1)) = f1.a*sqrt(t1)+f1.b;
end
hold on
h3 = plot(t2,z,':k','linewidth',2);
p_opd = p;


depths = 1.3;
t4text = max(max([times22,tGauss2]))+1;
text(t4text,depths(end)-0.4,['R^2 = ', num2str(R(3)), ' %'],'Color','black')
text(t4text,depths(end)-0.2,['R^2 = ', num2str(R(2)), ' %'],'Color','blue')
text(t4text,depths(end),['R^2 = ', num2str(R(1)), ' %'],'Color','red')
legend([h1 h2 h3],'TSR','NTC', 'OPD')%,'Location','northwest')
set(get(gca,'xlabel'),'string','Time max (s)','fontsize',12)
set(get(gca,'ylabel'),'string','Depth (mm)','fontsize',12)
set(get(gca,'title'),'string','Depth estimation','fontsize',11, 'FontWeight','bold')

depths = 0.2:0.2:1.8;
Zestimated;
Zexact = [depths;depths;depths];
Zerror = (100*abs(Zestimated-Zexact)./Zexact);
clc
disp('*********** Predicted values (mm) *******')
disp('Value     TSR         NTC         OPD')
disp(num2str([depths(1:6)' Zestimated(:,1:6)']))
disp('*********** Relative errors (%) *******')
disp('  TSR           NTC           OPD')
disp(num2str(Zerror(:,1:6)'))

    
   
   