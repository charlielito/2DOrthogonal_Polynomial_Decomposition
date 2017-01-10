close all; clear all; clc;
% Add to path 'functions' directory from whereever you run this 
% (obviosly that folder must be with this file)
currdir = mfilename('fullpath');
sep = strfind(currdir,'\');
path(path,[currdir(1:sep(end)),'functions']);


disp('-------- Loading data... ------')
load CFRPsample %RTfront tfront RTrear trear
disp('-------- Data loaded! ------')

[m,n,p1] = size(RTfront);
[m2,n2,p2] = size(RTrear);
p = min([p1 p2]);

% Front Sample defects Centers
mc = [73   152   228   304   380];
nc = [85   164   238   314   394];
%to fit better centers due to deformities of material
offsetm  = [0 -2 0 3 0; 2 -3 0 0 2; 2 0 0 0 2; 4 -2 2 0 4; 8 0 0 0 8];
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
kk = 2;
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

PCT = zeros(5,5,p1);
PCT2 = zeros(5,5,p2);

disp('-------- Calculating TSR ... (this may take some minutes) ------')
FTfront = TSR(RT_coefffront,tfront,9,'FT');
disp('-------- TSR front terminated! ... ------')
FTrear = TSR(RT_coeffrear,trear,9,'FT');
disp('-------- TSR terminated! ... ------')

disp('-------- Calculating OPD for biggest defects... ------')
h = waitbar(0,'Calculating OPD for biggest defects...');
steps = p;
coord=[5 1; 1 2; 5 3; 1 4; 5 5]; %only the biggest defects
for pp = 1:p
waitbar(pp / steps)
    for coor = coord'
        mm = coor(1);
        nn = coor(2);
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
close(h)
disp('-------- OPD calculated! ------')
%%
disp('-------- Calculating Normalized Contrast... ------')
NCfront = ntc(RTfront,tfront,5,100,380,mc,nc);
%NCfront = ntc(RTfront,tfront,5,184,74,mc,nc);
disp('-------- Normalized Contrast for front part terminated! ------')
NC = ntc(RTrear,trear,5,90,360,mc2,nc2); %for rear sample
%NC = ntc(RTrear,trear,5,54,360,mc2,nc2);
disp('-------- Normalized Contrast for REAR part terminated! ------')


%%
disp('-------- Calculating PCT for biggest defects... ------')
h = waitbar(0,'Calculating PCT for biggest defects...');
i = 0;
for coor = coord'
    mm = coor(1);
    nn = coor(2);
    waitbar(i / length(coord))
    intM = mi(mm,nn):mf(mm,nn);
    intN = ni(mm,nn):nf(mm,nn);
    intM2 = mi2(mm,nn):mf2(mm,nn);
    intN2 = ni2(mm,nn):nf2(mm,nn);

    Matrix_mn = squeeze(FTfront(intM,intN,:));
    Matrix_mn2 = squeeze(FTrear(intM2,intN2,:));
    PCA_Mat = compress3DTo2DForPCA(Matrix_mn);
    [~,~,V] = svd(PCA_Mat);
    PCT(mm,nn,:) = V(:,2)'; %time profile from 2nd orthogonal basis

    PCA_Mat = compress3DTo2DForPCA(Matrix_mn2);
    [~,~,V] = svd(PCA_Mat);
    PCT2(mm,nn,:)  =  V(:,2)'; %time profile from 2nd orthogonal basis
    i = i+1;
end
waitbar(i / length(coord))
close(h)
disp('-------- PCT calculated! ------')

%%
Tsfront = fitSoundA(RTfront,tfront,5,184,74); %sound areas fitted for PPT
Tsrear = fitSoundA(RTrear,trear,5,54,360);

Tsfront = fitSoundA(RTfront,tfront,5,184,74);
Tsrear = fitSoundA(RTrear,trear,5,54,360);

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
PCTFront = [];
PCTBack = [];
for nn = 1:5
    %OPD
    CoeffsFront = [CoeffsFront; squeeze(K0mn( Sm(mm,Hn(nn)), Hn(nn), :))'];
    CoeffsBack = [CoeffsBack; squeeze(K0mn2( Sm(mm,Hn2(nn)), Hn2(nn), :))'];
    
    %PCT
    PCTFront = [PCTFront; squeeze(PCT( Sm(mm,Hn(nn)), Hn(nn), :))'];
    PCTBack = [PCTBack; squeeze(PCT2( Sm(mm,Hn2(nn)), Hn2(nn), :))'];
    
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

Nn = min([length(PCTFront), length(PCTBack)]);
PCTTotal = [PCTFront(:,1:Nn); PCTBack(:,1:Nn)];

Nn = min([length(RTbigFront), length(RTbigBack)]);
RTbig = [RTbigFront(:,1:Nn); RTbigBack(:,1:Nn)];

tall = tfront(1:Nn);
% PPT Parameters
[~, frames] = size(RTbig);
Fs = 55;
nppt = 1;
Tsfront = Tsfront(1:nppt:frames);
Tsrear = Tsrear(1:nppt:frames);
Ts = [Tsfront;Tsfront;Tsfront;Tsfront;Tsfront;Tsrear;Tsrear;Tsrear;Tsrear;Tsrear];

depths = [0.2:0.2:1 1:0.2:1.8];

N = length(CoefTotal(1,:));
t5 = log((1:N)*fs);
Ndefect = 10;
tGauss = zeros(1,Ndefect);
tPol = zeros(1,Ndefect);

t = log((1:frames)*fs);
t3 = log(1:frames);
logt = log(tall);
times = zeros(1,10);
times2 = zeros(1,10);
times4 = zeros(1,10);
times5 = zeros(1,10);
for i = 1:Ndefect
    % *********** OPD calculations **********
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
    
    %************** NTC *******************
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
    t3 = log((2:frames)*fs);
    t3 = t3(1:end);
    p = polyfit(t3',Tdef',8);
    p = diffPoly(p);
    DTdef = polyval(p,t3);

    % *********** TSR **************************8
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
    xlim([min(log((1:frames)*fs)) max(t3) ])
    
    %************ PCT ************
    figure(4)
    subplot(2,5,i)
    C = PCTTotal(i,:);
    [Max, Idx] = findpeaks(C,'SORTSTR','ascend','minpeakdistance',100);
    DataInv = 1.01*max(C) - C;
    [Minima,MinIdx] = findpeaks(DataInv,'SORTSTR','ascend','minpeakdistance',100);
    plot(logt,C,'k')
    xlabel('ln(t)')
    ylabel('2nd Orthogonal Vector')
    hold on
    title(['Depth: ',num2str(depths(i)),' mm'])
    xlim([min(logt) max(logt) ])
    if isempty(Max) || Idx(end) < 15 || Idx(end) > 300% no maxima or fake maxima
    times4(i) = MinIdx(end);
    plot(logt(MinIdx(end)), C(MinIdx(end)),'rx')
    else
    times4(i) = Idx(end);
    plot(logt(Idx(end)),Max(end),'rx')
    end

    
    
    %********** PPT *****************
    figure(5)
    Tdef = RTbig(i,:);        
    y = log(Tdef);
    ypol = polynomialfit(logt(:),y(:),9);
    yfit = polyval(ypol,logt);
    Tdef2 = exp(yfit(1:nppt:frames));

    TdefFreq = fft(Tdef2);
    TsFreq = fft(Ts(i,:));
    Frames = size(TdefFreq,2);

    TdefFreq = TdefFreq(1:floor(Frames/2)+1);
    TsFreq = TsFreq(1:floor(Frames/2)+1);%only positive frequencies
    f = (0:Frames/2)*Fs/Frames/nppt;
    Pdef = angle(TdefFreq); %phase profiles
    Psound = angle(TsFreq);
    PCon = Psound-Pdef; %phase contrast
    factor = 10; %resampling factor to have better blind frequency estimation
    PCon2 = resample(double(PCon), factor,1);
    f2 = (0:length(PCon2)-1)*Fs/Frames/nppt/factor; %frequency resampled
    times5(i) = f2(zerocrossPPT(PCon2(2:end))); %calculate zero cross (blind f_b)
    plot(f,PCon)
    xlabel('Frequency (Hz)')
    ylabel('Phase Contrast (rad)')
    hold all
    title('Phase contrasts for different dephts ')
    xlim([min(f) max(f) ])
end
figure(5),legend('0.2 mm', '0.4 mm', '0.6 mm','0.8 mm','1.0 mm','1.0 mm','1.2 mm','1.4 mm','1.6 mm','1.8 mm')

%% 
% ***** TSR *********************
depths = 0.2:0.2:1.8;
limit = 4;
t1 = times(1:limit);
figure(6), plot(t1, depths(1:limit),'xr', 'linewidth',1.5,'MarkerSize',10), hold on

degree = 2;
tf = 6;
[p,s] = polyfit(t1', depths(1:limit)',degree);
t2 = linspace(t1(1),tf);
z = polyval(p,t2);
h1 = plot(t2,z,'r','linewidth',1.5);
R = zeros(1,5);
R(1) =  (1 - s.normr^2/norm(depths(1:limit)-mean(depths(1:limit)))^2)*100;
Zestimated = zeros(5,length(depths));
Zestimated(1,1:length(t1)) = polyval(p,t1);
p_tsr = p;


% *************** NTC ************
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


% ********* OPD ****************
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

% ********* PCT **********
limit = 6;
depths = 0.2:0.2:1.8;
depths = depths(1:limit);
linear = 1; %linear fit or not linear (sqrt(x))
times4t = [times4(1:4), times4(6:7)];
t1 = times4t(1:limit)*fs;
plot(t1, depths,'*g', 'linewidth',1.5,'MarkerSize',7)
t2 = linspace(0,5);
if (linear)
    [p,s] = polyfit(t1', depths',1);
    z = polyval(p,t2);
    R(4) =  (1 - s.normr^2/norm(depths-mean(depths))^2)*100;
    Zestimated(4,1:length(t1)) = polyval(p,t1);
    p_pca = p;
else
    Eqn = 'a*sqrt(x) + b';
    [f2, gof] = fit(t1',depths',Eqn,'Startpoint',[1 1]);
    z = f2.a*sqrt(t2)+f2.b;
    R(4) = gof.rsquare*100;
    Zestimated(4,1:length(t1)) = f2.a*sqrt(t1)+f2.b;
    p_pca(1) = f2.a;
    p_pca(2) = f2.b;
end
hold on
h4 = plot(t2,z,'--g','linewidth',1.5,'markersize',3);
depths = 0.2:0.2:1.8;


% **** PPT *************
times5t = [times5(1:5), times5(7:end)];
range = 1:6;
t = times5t(range).^(-1);
plot(t, depths(range),'^m', 'linewidth',1.5,'MarkerSize',7)
[f1,gof] = fit(t',depths(range)','a*sqrt(x) + b');
t2 = linspace(0,6);
z = f1.a*sqrt(t2) + f1.b;
h5 = plot(t2,z,'m');
R(5) =  gof.rsquare*100;
p_ppt = [f1.a f1.b];
Zestimated(5,1:length(times5t)) = p_ppt(1)*sqrt(times5t.^-1) + p_ppt(2);


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

depths = 1.3;
t4text = max(max([times22,tGauss2]))+1.5;

text(t4text,depths(end)-0.42,['R^2_{TSR} = ', num2str(R(1)), ' %'],'Color','red','FontSize',8)
text(t4text,depths(end)-0.48,['z_{TSR} = ', num2str(p_tsr(1)), 't^2 ', ...
    signn(sign(p_tsr(2))),' ', num2str(abs(p_tsr(2))), 't ', ...
    signn(sign(p_tsr(3))),' ',  num2str(abs(p_tsr(3))) ],'Color','red','FontSize',8)

text(t4text,depths(end)-0.58,['R^2_{NTC} = ', num2str(R(2)), ' %'],'Color','blue','FontSize',8)
text(t4text,depths(end)-0.64,['z_{NTC} = ', num2str(p_ntc(1)), 't ', ...
    signn(sign(p_ntc(2))),' ',  num2str(abs(p_ntc(2))) ],'Color','blue','FontSize',8)

text(t4text,depths(end)-0.74,['R^2_{OPD} = ', num2str(R(3)), ' %'],'Color','black','FontSize',8)
text(t4text,depths(end)-0.8,['z_{OPD} = ', num2str(p_opd(1)), 't ', ...
    signn(sign(p_opd(2))),' ',  num2str(abs(p_opd(2))) ],'Color','black','FontSize',8)

text(t4text,depths(end)-0.9,['R^2_{PCT} = ', num2str(R(4)), ' %'],'Color','g','FontSize',8)
text(t4text,depths(end)-.96,['z_{PCT} = ', num2str(p_pca(1)), 't ', ...
    signn(sign(p_pca(2))),' ',  num2str(abs(p_pca(2))) ],'Color','g','FontSize',8)

text(t4text-1.7,depths(end)-0.9,['R^2_{PPT} = ', num2str(R(5)), ' %'],'Color','m','FontSize',8)
text(t4text-1.7,depths(end)-0.96,['z_{PPT} = ', num2str(p_ppt(1)), 't^{1/2} ', ...
    signn(sign(p_ppt(2))),' ',  num2str(abs(p_ppt(2))) ],'Color','m','FontSize',8)

legend([h1 h2 h3 h4 h5],'TSR','NTC', 'OPD','PCT', 'PPT','Location','northwest')
set(get(gca,'xlabel'),'string','Specific characteristic time t (s)','fontsize',12)
set(get(gca,'ylabel'),'string','Depth z (mm)','fontsize',12)
set(get(gca,'title'),'string','Depth estimation','fontsize',11, 'FontWeight','bold')
axis( [0 5.7 0.1 1.4])


depths = 0.2:0.2:1.8;
Zestimated;
Zexact = [depths;depths;depths;depths;depths];
Zerror = (100*(round(100*Zestimated)/100-Zexact)./Zexact);
Zerror = Zerror(:,1:6);
clc
disp('*********** Predicted values (mm) *******')
disp('Value        TSR        NTC         OPD         PCT         PPT')
disp(num2str([depths(1:6)' Zestimated(:,1:6)'],'%0.2f        '))
disp('*********** Relative errors (%) *********')
disp('            TSR           NTC         OPD          PCT          PPT')
disp(num2str([nan*ones(1,size(Zerror,2))' Zerror'],'%0.1f       '))