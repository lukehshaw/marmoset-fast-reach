
function pursuit_model_compare_fpm()

% Generates Figure 3G, Figure 3H from  
% Shaw,L, Wang KH, Mitchell, J (2024) Fast Prediction in Marmoset Reach-to-Grasp Movements for Dynamic Prey.
%
% Population range vector derivative corrlation to range vector using reaching data and pure pursuit and proportional navigation
% simulations (3G). Target-congruent hand lateral velocity population average for reaches for reaching data and simulations.
% Luke Shaw, Kuan Hong Wang, and Jude Mitchell 4/2023
% MATLAB R2022b
%
% Reaching data structure marmo_reach_model.mat available at
% https://doi.org/10.5281/zenodo.7869286
%
% Range vector derivative correlation analysis modified from:
% Mischiati, M., et al. (2015). "Internal models direct dragonfly interception steering." Nature 517(7534): 333-338.
% Pursuit simulation scripts modified from:
% Brighton, C. H. and G. K. Taylor (2019). "Hawks steer attacks using a guidance system tuned for close pursuit of erratically manoeuvring targets." Nat Commun 10(1): 2462.

%%

load marmo_reach_model.mat
nReach = size(model.x.cricket,2);
vecCorA = cell(nReach,1);
diffvecCorA = cell(nReach,1);
diffvecCorB = cell(nReach,1);
%***
velvecCorA = cell(1,3);
for kk = 1:3
    velvecCorA{kk} = cell(nReach,1);
end
%*****
nq = 31; %length of interpolated reaches
smw = 10; % smooth window

downSample=1;

x1=[0,0]; %create arena outline
x2=[9.65,8.55];
x3=[19.2,0];
x4=[2.8,5.3];
x5=[16.2,5.53];
Arena=polyfit([x1(1) x2(1) x3(1) x4(1) x5(1)],[x1(2) x2(2) x3(2) x4(2) x5(2)],4);

bestSims = [];
lamb = 1.0;  % 0 - weight position fit, 1 weight velocity fit 
             % using position error, R2 is PP: 89.1(1.5)
             %                             PN: 93.3(1.3)
             %                             PNP: 93.6(1.1)
             %                             PND: 90.2(1.9) - 50 ms delay
             % using velocity error, R2 is PP: 95.0(0.6)
             %                             PN: 96.0(0.5) ... 
minoff = 18; % force identical following this number of steps
DTau = 18; % consider PP and PN models forced at tau=12 is 50ms
shouldVis=0; %1=yea vis 0=no vis

%******* fit entire dataset with one parameterization
bestPP_R2 = 0;
bestPP_gain = NaN;
%*******
bestPN_R2 = 0;
bestPN_gain = NaN;
%******
bestPND_R2 = 0;
bestPND_gain = NaN;
%*********
Gmin = 0.1;
Gmax = 50;
Gn = 20;
GSAMP = exp( log(Gmin):((log(Gmax)-log(Gmin))/Gn):log(Gmax) );
for gain = GSAMP % Brighton + Taylor, K or N term
    disp(['Stepping ',num2str(gain)]);  
    if (1) % used to be loop over Tau
         % disp(['Tau: ',num2str(tau)]);    
         %********************        
         PPsum = 0; sumPPn = 0;
         PNsum = 0; sumPNn = 0;
         PNDsum = 0; sumPNDn = 0;
         for i = 1:nReach % number of reaches
           Hx= smoothdata(model.x.hand{i}(1:end,:),'gaussian',smw);
           Hy= smoothdata(model.y.hand{i}(1:end,:),'gaussian',smw);
           Cx = vertcat(model.x.cricketexfull{i},model.x.cricket{i});
           Cy = vertcat(model.y.cricketexfull{i},model.y.cricket{i});
           Cx = smoothdata(Cx,'gaussian',smw);
           Cy = smoothdata(Cy,'gaussian',smw);
           
           %******allow trial by trial start direction search (to account
           %****** for arbitrary start posture with each trial)
           [PPx,PPy,PPr2,PPn] = FindBestPNP(Hx,Hy,Cx,Cy,gain,0,1,lamb,minoff,0);
           [PNx,PNy,PNr2,PNn] = FindBestPNP(Hx,Hy,Cx,Cy,0,gain,1,lamb,minoff,0);       
           [PNDx,PNDy,PNDr2,PNDn] = FindBestPNP(Hx,Hy,Cx,Cy,0,gain,DTau,lamb,minoff,0);       
           %****
           if ~isnan(PPr2)
             PPsum = PPsum + (PPn * PPr2);
             sumPPn = sumPPn + PPn;
           end
           if ~isnan(PNr2)
             PNsum = PNsum + (PNn * PNr2);
             sumPNn = sumPNn + PNn;
           end
           if ~isnan(PNDr2)
             PNDsum = PNDsum + (PNDn * PNDr2);
             sumPNDn = sumPNDn + PNDn;
           end
         end
         PPsum = (PPsum/sumPPn);  % norm over trials, total R2
         PNsum = (PNsum/sumPNn);  % norm over trials, total R2
         PNDsum = (PNDsum/sumPNDn);  % norm over trials, total R2
         
         %*********************
         if (PPsum > bestPP_R2)
              bestPP_gain = gain;
              bestPP_R2 = PPsum;
         end
         %*****************
         if (PNsum > bestPN_R2)
              bestPN_gain = gain;
              bestPN_R2 = PNsum;
         end
         %*****************
         if (PNDsum > bestPND_R2)
              bestPND_gain = gain;
              bestPND_R2 = PNDsum;
         end
         %*********
   end
end

%**********

%********* fit a mixed model, but with a delay
bestPF_R2 = 0;
bestPF_K = NaN;
bestPF_N = NaN;
%******
for K = GSAMP % max 200 ms into future
   disp(['Stepping K ',num2str(K)]);  
   for N = GSAMP 
         %********************        
         PFsum = 0;
         sumPFn = 0;
         for i = 1:nReach % number of reaches           
           %***** NOTE you need to smooth with them concatenated
           Hx= smoothdata(model.x.hand{i}(1:end,:),'gaussian',smw);
           Hy= smoothdata(model.y.hand{i}(1:end,:),'gaussian',smw);
           Cx = vertcat(model.x.cricketexfull{i},model.x.cricket{i});
           Cy = vertcat(model.y.cricketexfull{i},model.y.cricket{i});
           Cx = smoothdata(Cx,'gaussian',smw);
           Cy = smoothdata(Cy,'gaussian',smw);
           
           %********
           [PFx,PFy,PFr2,PFn] = FindBestPNP(Hx,Hy,Cx,Cy,K,N,1,lamb,minoff,0);
           if ~isnan(PFr2)
             PFsum = PFsum + (PFn*PFr2);
             sumPFn = sumPFn + PFn;
           end
         end
         PFsum = PFsum/sumPFn;
         
         %******************
         if (PFsum > bestPF_R2)
              bestPF_K = K;
              bestPF_N = N;
              bestPF_R2 = PFsum;
         end
         %*********
   end
end
%*************************          


for i = 1:nReach % number of reaches

    nFrame = length(model.x.hand{i});
    P = linspace(1,nFrame,nq);
    x = (1:nFrame)';
    xv = (1:nFrame-1)';

    %***** NOTE you need to smooth with them concatenated
    zCx = vertcat(model.x.cricketexfull{i},model.x.cricket{i});
    zCy = vertcat(model.y.cricketexfull{i},model.y.cricket{i});
    zCx = smoothdata(zCx,'gaussian',smw);
    zCy = smoothdata(zCy,'gaussian',smw);
    %*******       
    Hx= model.x.hand{i}(1:end,:);
    Hy= model.y.hand{i}(1:end,:);
    Cx= model.x.cricket{i}(1:end,:);
    Cy= model.y.cricket{i}(1:end,:);
    Cx2 = model.x.cricketexfull{i};
    Cy2 = model.y.cricketexfull{i};
    Hx = smoothdata(Hx,'gaussian',smw);
    Hy = smoothdata(Hy,'gaussian',smw);
    M1 = length(Cx);
    M2 = length(Cx2);
    Sx = smoothdata(vertcat(Cx2,Cx),'gaussian',smw);
    Sy = smoothdata(vertcat(Cy2,Cy),'gaussian',smw);
    Cx2 = flipud(Sx(1:M2));
    Cy2 = flipud(Sy(1:M2));
    Cx = Sx((M2+1):(M2+M1));
    Cy = Sy((M2+1):(M2+M1));

    
    [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1] = kinematicAnalysisFutureZ_edit([Hx Hy],[Cx Cy],1);

    %General Velocity
    veloHandtotal{i}=vecnorm([VHo1,VHp1],2,2).*240;
    % veloHandortho{i}=abs(VHo1).*240;
    veloHandortho{i} = (VHo1 .* sign(VCo1)) .* 240;  % set by cricket sign
    veloCrickettotal{i}=vecnorm([VCo1,VCp1],2,2).*240;
    veloCricketortho{i}=abs(VCo1).*240;
        
    rangeVectors = [Cx-Hx,Cy-Hy];
    % nFrame x 2 matrix, columns are x and y coordinates
    % derivative of range vectors
    dRV = diff(rangeVectors);
    dRVjm=veloHandtotal{i}-veloCrickettotal{i};

        
    %******* plot the trajectory of trial with range vectors
    %****** find best pure-pursuit to match observed trajectory
    %****** given the same hand speed moment to moment (only
    %****** differing by the steering
    %[PPx,PPy,PPerr] = FindBestPurePursuit(Hx,Hy,Cx,Cy,...
    %                      bestPP_gain,bestPP_tau,lamb,minoff);
    [PPx,PPy,PPr2,n] = FindBestPNP(Hx,Hy,zCx,zCy,...
                          bestPP_gain,0,1,lamb,minoff,0);
    [PNx,PNy,PNr2,n] = FindBestPNP(Hx,Hy,zCx,zCy,...
                          0,bestPN_gain,1,lamb,minoff,0);
    [PFx,PFy,PFr2,n] = FindBestPNP(Hx,Hy,zCx,zCy,...
                          bestPF_K,bestPF_N,1,lamb,minoff,0);
    %** PN but forced to same delay (good on straight, but unstable
    %** when the turns get larger
    [PNDx,PNDy,PNDr2,n] = FindBestPNP(Hx,Hy,zCx,zCy,...
                          0,bestPND_gain,DTau,lamb,minoff,0); 
    %***********
    
    %******* velocity for hand models
    [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1] = kinematicAnalysisFutureZ_edit([PPx PPy],[Cx Cy],1);
    PPveloHandtotal{i}=vecnorm([VHo1,VHp1],2,2).*240;
    PPveloHandortho{i} = (VHo1 .* sign(VCo1)) .* 240;  % set by cricket sig
    %******* computer lateral velocities for models
    [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1] = kinematicAnalysisFutureZ_edit([PNx PNy],[Cx Cy],1);
    PNveloHandtotal{i}=vecnorm([VHo1,VHp1],2,2).*240;
    PNveloHandortho{i} = (VHo1 .* sign(VCo1)) .* 240;  % set by cricket sig
    %***********
    
    if (shouldVis == 1)
            hf = figure(2); hold off;
            set(hf,'Position',[100 200 900 500]);
            skip = 4;
            rcolo = [1,0.5,0.5];
            bcolo = [0,0.5,0.5];
            gcolo = [0.7,0.7,0.7];
            for tt = (1+skip):skip:size(rangeVectors,1)
              plot([Hx(tt),Cx(tt)],[Hy(tt),Cy(tt)],'k-','Color',gcolo); hold on;
              plot([Hx(tt),Hx(tt-skip)],[Hy(tt),Hy(tt-skip)],'b-','Color',bcolo);
              plot(Hx(tt),Hy(tt),'b.','Color',bcolo,'Markersize',8); 
              plot([Cx(tt),Cx(tt-skip)],[Cy(tt),Cy(tt-skip)],'r-','Color',rcolo);
              plot(Cx(tt),Cy(tt),'r.','Color',rcolo,'Markersize',8); 
              %******
              plot([PPx(tt),PPx(tt-skip)],[PPy(tt),PPy(tt-skip)],'k-');
              %**********
              plot([PNx(tt),PNx(tt-skip)],[PNy(tt),PNy(tt-skip)],'b-');
              %**********
              plot([PNDx(tt),PNDx(tt-skip)],[PNDy(tt),PNDy(tt-skip)],'b:');
              %**********
              plot([PFx(tt),PFx(tt-skip)],[PFy(tt),PFy(tt-skip)],'k-','Color',[0.8,0.4,0]);
              %*******
            end
            axis([0 20 -2 10]);
            ht = title(sprintf(['T(%d) M(R2,gain,tau) PP',...
                           '(%4.3f,%3.1f,%d) PN(%4.3f,%3.1f,%d)',...
                           ' PF(%4.3f,%3.1f,%3.1f,%d) PND(%4.3f,%3.1f,%d)'],...
                      i,PPr2,bestPP_gain,1,...
                        PNr2,bestPN_gain,1,...
                        PFr2,bestPF_K,bestPF_N,1,...
                        PNDr2,bestPND_gain,floor(DTau)));
            set(ht,'Fontsize',10);
            [PPr2,PNr2,PFr2,PNDr2]
            input('check');
            %***********
        end
        bestSims = [bestSims ; [PPr2,PNr2,PFr2,PNDr2] ];
        
        % correlation between rangeVectors and their derivative vectors
        RV = rangeVectors(1:end-1,:); % to match with derivative row number
        RVnorm = vecnorm(RV,2,2);
        dRVnorm = vecnorm(dRV,2,2);
        vecCor = sum(RV.*dRV,2)./(RVnorm.*dRVnorm);
        vecCorA{i} = vecCor;

        %******** same correlation, but for the fit pure pursuit trajectories
        if (1)
           %***** range * d range, pure pursuit 
           PPrangeVectors = [Cx-PPx,Cy-PPy];  % get range vec of pure pursuit
           dPPrangeVectors = diff(PPrangeVectors);
           RV = PPrangeVectors(1:end-1,:); % to match with derivative row number
           RVnorm = vecnorm(RV,2,2);
           dRVnorm = vecnorm(dPPrangeVectors,2,2);
           diffvecCor = sum(RV.*dRV,2)./(RVnorm.*dRVnorm);
           diffvecCorA{i} = diffvecCor;

           %**** range * d range, parallel navigation
           PPrangeVectors = [Cx-PNx,Cy-PNy];  % get range vec of pure pursuit
           dPPrangeVectors = diff(PPrangeVectors);
           RV = PPrangeVectors(1:end-1,:); % to match with derivative row number
           RVnorm = vecnorm(RV,2,2);
           dRVnorm = vecnorm(dPPrangeVectors,2,2);
           diffvecCor = sum(RV.*dRV,2)./(RVnorm.*dRVnorm);
           diffvecCorB{i} = diffvecCor;
        end
        
        %******* get correlations with actual to model velocities
        AN = 1;
        RV = [diff(Hx,AN),diff(Hy,AN)];  % hand velocity
        RVlat = zeros(size(RV)); % 1d vector
        %*** project onto range vectors to get lateral component
        if (0)
          R = [(Cx-Hx),(Cy-Hy)];
          for tt = 1:size(RV,1)
            RR = R(tt,:)/norm(R(tt,:));
            RVlat(tt,1) = RV(tt,1)*RR(2) - RV(tt,2)*RR(1);  % lateral comp
          end
          RV = RVlat;
        end
        %********
        RVnorm = vecnorm(RV,2,2);  
        for kk = 1:3
          if (kk == 1)
             dRV = [diff(PPx,AN),diff(PPy,AN)];
          end
          if (kk == 2)
             dRV = [diff(PNx,AN),diff(PNy,AN)];
          end
          if (kk == 3)
             dRV = [diff(PFx,AN),diff(PFy,AN)];
          end
          %********
          %*** project onto range vectors to get lateral component
          if (0)
            R = [(Cx-Hx),(Cy-Hy)];
            dRVlat = zeros(size(RV)); % 1d vector
            for tt = 1:size(dRV,1)
              RR = R(tt,:)/norm(R(tt,:));
              dRVlat(tt,1) = dRV(tt,1)*RR(2) - dRV(tt,2)*RR(1);  % lateral comp
            end
            dRV = dRVlat;
          end
          %*******
          dRVnorm = vecnorm(dRV,2,2);
          velvecCor = sum(RV.*dRV,2)./(RVnorm.*dRVnorm);
          velvecCorA{kk}{i} = velvecCor;
        end
end

bestSims
nanmean(bestSims)
nanstd(bestSims)/sqrt(sum(~isnan(bestSims(:,1))))
zbestSims = bestSims(:,2:end);
%*********
zz = find(~isnan(bestSims(:,1)));
p1 = signrank(bestSims(zz,1),bestSims(zz,2));  % PP vs PN model
p2 = signrank(bestSims(zz,2),bestSims(zz,3));  % PN vs PNP model
p3 = signrank(bestSims(zz,2),bestSims(zz,4));  % PN vs PNP model

%*** parameters ********
disp('Sign rank test (PP vs PN), (PN vs PNP), (PN vs PNdelay) based on trial by trial R2');
[p1,p2,p3]
disp(sprintf('M(R2,gain,tau) PP(%3.2f,%d) PN(%3.2f,%d)',...
              bestPP_gain,1,bestPN_gain,1));
disp(sprintf('               PF(%3.2f,%3.2f,%d) PND(%3.2f,%d)',...
              bestPF_K,bestPF_N,1,bestPND_gain,floor(DTau)));
%******            

for k = 1:size(zbestSims,1)
    for kk = 1:size(zbestSims,2)
        zbestSims(k,kk) = 100 * (zbestSims(k,kk) - bestSims(k,1)) / ...
                                (zbestSims(k,kk) + bestSims(k,1));
    end        
end
[0, nanmean( zbestSims )]
input('see all');

%% plot final results
%******* summary stats on reaching movements (Luke's original code)
veloHand=veloHandortho;
veloCricket=veloCricketortho;
vcSize=max(cellfun('size',veloCricket,1)); %find cell with longest length
vcAC=NaN(size(veloCricket,2),vcSize);
for q=1:size(veloCricket,2)
    vcAC(q,(vcSize+1)-length(veloCricket{q}):end)=veloCricket{q}';
end
vcSize=max(cellfun('size',veloHand,1)); %find cell with longest length
vcAH=NaN(size(veloHand,2),vcSize);
for q=1:size(veloHand,2)
        vcAH(q,(vcSize+1)-length(veloHand{q}):end)=veloHand{q}';
end

%******* get lateral velocity for best PP model ************
PPveloHand = PPveloHandortho;
PPvcSize=max(cellfun('size',PPveloHand,1)); %find cell with longest length
PPvcAH=NaN(size(PPveloHand,2),PPvcSize);
for q=1:size(PPveloHand,2)
        PPvcAH(q,(PPvcSize+1)-length(PPveloHand{q}):end)=PPveloHand{q}';
end
%****
PNveloHand = PNveloHandortho;
PNvcSize=max(cellfun('size',PNveloHand,1)); %find cell with longest length
PNvcAH=NaN(size(PPveloHand,2),PNvcSize);
for q=1:size(PNveloHand,2)
        PNvcAH(q,(PNvcSize+1)-length(PNveloHand{q}):end)=PNveloHand{q}';
end
%********


%***
figure
for i = 1:nReach
    tmp = vecCorA{i};
    nFrame = length(tmp);
    plot((-nFrame+1:0)/.240,vecCorA{i},'Color',[.8 .8 .8]);
    hold on;
end

vcSize=max(cellfun('size',vecCorA,1)); %find cell with longest length
vcA=NaN(size(vecCorA,1),vcSize);
for q=1:size(vecCorA,1)
    vcA(q,(vcSize+1)-length(vecCorA{q}):end)=vecCorA{q}';
end

vcSize2=max(cellfun('size',diffvecCorA,1)); %find cell with longest length
diffvcA=NaN(size(diffvecCorA,1),vcSize2);
for q=1:size(diffvecCorA,1)
    diffvcA(q,(vcSize2+1)-length(diffvecCorA{q}):end)=diffvecCorA{q}';
end

vcSize3=max(cellfun('size',diffvecCorB,1)); %find cell with longest length
diffvcB=NaN(size(diffvecCorB,1),vcSize3);
for q=1:size(diffvecCorB,1)
    diffvcB(q,(vcSize3+1)-length(diffvecCorB{q}):end)=diffvecCorB{q}';
end

velvcA = cell(1,3);
for kk = 1:3
    velSize(kk)=max(cellfun('size',velvecCorA{kk},1)); %find cell with longest length
    velvcA{kk}=NaN(size(velvecCorA{kk},1),velSize(kk));
    for q=1:size(velvecCorA{kk},1)
       velvcA{kk}(q,(velSize(kk)+1)-length(velvecCorA{kk}{q}):end)=velvecCorA{kk}{q}';
    end
end

xx = (-vcSize+1:0)/.240;
uu = nanmean(vcA);
normo = sum(~isnan(vcA));
su = nanstd(vcA) ./ sqrt(normo);
%**************
%plot((-vcSize+1:0)/.240,nanmean(vcA,1),'k-','linewidth',4,'Color',[1,0,1]);
aa = [xx fliplr(xx)];
bb = [(uu-su) fliplr(uu+su)];
fill(aa,bb,[1,0,1],'FaceAlpha',0.2,'Linestyle','none'); hold on;
plot(xx,uu,'k-','linewidth',4,'Color',[0.8,0,0.8]);

if (1)
  xx = (-vcSize2+1:0)/.240;
  uu2 = nanmean(diffvcA);
  normo2 = sum(~isnan(diffvcA));
  su2 = nanstd(diffvcA) ./ sqrt(normo2);
  aa = [xx fliplr(xx)];
  bb = [(uu2-su2) fliplr(uu2+su2)];
  %fill(aa,bb,[1,0,0],'FaceAlpha',0.2,'Linestyle','none');
  plot(xx,uu2,'k--','linewidth',2,'Color',[0,0,0]);

  xx = (-vcSize2+1:0)/.240;
  uu2 = nanmean(diffvcB);
  normo2 = sum(~isnan(diffvcB));
  su2 = nanstd(diffvcB) ./ sqrt(normo2);
  aa = [xx fliplr(xx)];
  bb = [(uu2-su2) fliplr(uu2+su2)];
  %fill(aa,bb,[1,0,0],'FaceAlpha',0.2,'Linestyle','none');
  plot(xx,uu2,'k--','linewidth',2,'Color',[0,0,1]);
end
        
ylim([-1 0])
xlim([-220 0])
set(gca,'TickLength',[0 0])
title('PP Real Time Series');
ylabel('Correlation');
xlabel('ms')

%********* new figure to plot correlations with models
figure;
colo = [[0,0,0];[0,0,1];[0.6,0.6,0]];
for kk = 1:3
  xx = (-velSize(kk)+1:0)/.240;
  uu2 = nanmean(velvcA{kk});
  normo2 = sum(~isnan(velvcA{kk}));
  su2 = nanstd(velvcA{kk}) ./ sqrt(normo2);
  aa = [xx fliplr(xx)];
  bb = [(uu2-su2) fliplr(uu2+su2)];
  fill(aa,bb,colo(kk,:),'FaceAlpha',0.2,'Linestyle','none'); hold on;
  plot(xx,uu2,'k--','linewidth',2,'Color',colo(kk,:));
end        
ylim([0.5 1])
xlim([-220 0])
set(gca,'TickLength',[0 0])
title('Time Series with Models');
ylabel('Correlation');
xlabel('ms')

%******* plot the mean velocities (direct and lateral) over reaches
figure;
xx = ((-vcSize+2):0)/.240;
plot(xx,zeros(size(xx)),'k--'); hold on;
uuc = nanmean(vcAC,1);
suc = nanstd(vcAC,1) ./ sqrt( sum(~isnan(vcAC)) );
aa = [xx fliplr(xx)];
bb = [(uuc+suc) fliplr(uuc-suc)];
fill(aa,bb,[1,0.5,0.5],'FaceAlpha',0.3,'Linestyle','none');
plot(xx,uuc,'r-','Color',[1,0.5,0.5],'linewidth',2); hold on
%****
uuh = nanmean(vcAH,1);
suh = nanstd(vcAH,1) ./ sqrt( sum(~isnan(vcAH)) );
aa = [xx fliplr(xx)];
bb = [(uuh+suh) fliplr(uuh-suh)];
fill(aa,bb,[0,0.6,0.6],'FaceAlpha',0.3,'Linestyle','none');
plot(xx,uuh,'b-','Color',[0,0.6,0.6],'linewidth',2);

if (1)
%***** hand lateral speed for PP and PN models
  PPuuh = nanmean(PPvcAH,1);
  PPsuh = nanstd(PPvcAH,1) ./ sqrt( sum(~isnan(PPvcAH)) );
  aa = [xx fliplr(xx)];
  bb = [(PPuuh+PPsuh) fliplr(PPuuh-PPsuh)];
  % fill(aa,bb,[0,0,0],'FaceAlpha',0.3,'Linestyle','none');
  plot(xx,PPuuh,'b--','Color',[0,0,0],'linewidth',2);
  %***** hand lateral speed for PP model
  PNuuh = nanmean(PNvcAH,1);
  PNsuh = nanstd(PNvcAH,1) ./ sqrt( sum(~isnan(PNvcAH)) );
  aa = [xx fliplr(xx)];
  bb = [(PNuuh+PNsuh) fliplr(PNuuh-PNsuh)];
  % fill(aa,bb,[0,0,1],'FaceAlpha',0.3,'Linestyle','none');
  plot(xx,PNuuh,'b--','Color',[0,0,1],'linewidth',2);
end

%********
xlim([-220 0])
ylim([-5 15])
set(gca,'TickLength',[0 0])
title('Mean Velos Real Time Series');
ylabel('cm/s');
xlabel('ms')

return;


%********* BELOW, TRYING TO INCORPORATE USING EXTENDED CRICK VELOCITIES 
%Copyright. Graham K. Taylor & Caroline H. Brighton, University of Oxford, November 2018.
function [PPx,PPy,r2,N] = FindBestPNP(Hx,Hy,Cx,Cy,K,N,tau,lamb,minoff,taup)
      
      stepSize = 240;
      % tauMax = max(tau,minoff);
      tauMax = minoff;
      if (tauMax > (length(Hx)-1) )
          PPx = Hx;
          PPy = Hy;
          r2 = NaN;
          N = 0;
          return;
      end      
      birdPos = [Hx,Hy,zeros(size(Hx))];
      birdVel = [diff(Hx),diff(Hy),zeros(size(diff(Hx)))] * stepSize;
      birdVel = [birdVel(1,:) ; birdVel];
      birdSpeed = vecnorm(birdVel')';
      lurePos = [Cx,Cy,zeros(size(Cx))];
      lureVel = [diff(Cx),diff(Cy),zeros(size(diff(Cx)))] * stepSize;
      lureVel = [lureVel(1,:) ; lureVel];
      %****** call their routine
      simPos = zfitPNP(birdPos,birdVel,birdSpeed,lurePos,lureVel,stepSize,tauMax,K,N,tau,taup);
      %*******
      PPx = simPos(:,1);
      PPy = simPos(:,2);
      ex = tauMax:length(Hx);
      Ptot = sum( (Hx(ex)-mean(Hx(ex))) .^ 2) + sum( (Hy(ex)-mean(Hy(ex))) .^ 2 );
      Perr = sum( (Hx(ex) - PPx(ex)) .^ 2) + sum( (Hy(ex) - PPy(ex)) .^ 2);
      Vtot = sum( diff(Hx(ex)) .^ 2) + sum( diff(Hy(ex)) .^ 2);
      Verr = sum( (diff(Hx(ex)) - diff(PPx(ex))) .^ 2) + sum( (diff(Hy(ex)) - diff(PPy(ex))) .^ 2);    
      %*******
      N = (length(Hx)-tauMax);
      r2 = ((1-lamb)*((Ptot-Perr)/Ptot) + lamb * ((Vtot-Verr)/Vtot));   % weight position and vel R2
      %********
return


%Copyright. Graham K. Taylor & Caroline H. Brighton, University of Oxford, November 2018.
function [simPos] = zfitPNP(birdPos,birdVel,birdSpeed,lurePos,lureVel,stepSize,tauMax,K,N,tau,taup)

  %Initializes matrices for bird's simulated position and velocity
  simPos=birdPos*nan; simPos(1:tauMax,:)=birdPos(1:tauMax,:);
  simVel=birdVel*nan; % simVel(1:tauMax,:)=birdVel(1:tauMax,:);

  dN = size(lurePos,1) - size(birdPos,1);  % cricket trace is appended (sub dN to sync times)
  
  %*** rest to straight line for history up to TauMax (reduce noise)
  uvel = nanmean(birdVel(1:tauMax,:));
  for t = 1:tauMax
      simVel(t,:) = uvel;  % fix at least a few steps early to fixed hand vel (noise reduction)
  end
  %****
  cvel = nanmean(lureVel(1:tauMax,:));  % if negative time missing, use extrapolation back in time lure vel
  %*********
  
  for t=tauMax+1:length(birdPos)
    
    %% Uses PP guidance law to compute bird's commanded acceleration at time step (t-1)
      
        %Calculates line-of-sight vector (r) lagged by tau
        % r=lurePos(t-tau,:)-simPos(t-tau,:);
     
        %**** implement reverse approximation if cricket position missing
        if (0) % ((t-tau+dN) < 1)
            db = t-tau+dN-1;
            crickpos = lurePos(1,:) + ((db/stepSize) * cvel(1,:));
            crickvel = cvel(1,:);
        else
            crickpos = lurePos(t-tau+dN,:);
            crickvel = lureVel(t-tau+dN,:);
        end
        
        % implement reverse approximation on hand position if missing
        % negative times
        if (0) % ( (t-tau)<1)
            db = t-tau-1;
            handpos = simPos(1,:) + ((db/stepSize) * simVel(1,:));
            handvel = simVel(1,:);
        else
            handpos = simPos(t-tau,:);
            handvel = simVel(t-tau,:);
        end
        
        %******* prediction from cricket position
        lurePred = (taup * crickvel) + crickpos;     
        r = lurePred - handpos;  %effectively using hand position a delay
         
        %Defines velocity vector (v) lagged by tau
        v=handvel;  % simVel(t-tau,:);
    
        if (0)
          %Calculates deviation angle (delta) lagged by tau
          delta=acos(dot(r,v)/(norm(r)*norm(v)));
    
          %Creates unit vector normal to line of sight (r) and velocity (v)
          n=cross(r,v)/(norm(r)*norm(v)); n=n/norm(n); 
    
          %Calculates commanded acceleration under delayed PP guidance law
          birdAcc=-K*delta*cross(n,simVel(t-1,:));
      
        else
            
          %**** translate into 2D pursuit vector (forget steering)
          birdAcc = K * r * stepSize;
              
        end
        
        %% Uses PN guidance law to compute bird's commanded acceleration at time (t-1)
        %Calculates line-of-sight vector (r) lagged by tau

        %Calculates line-of-sight rate (lambdaDot) at time (t-delay)
        % lambdaDot=cross(r, (lureVel(t-tau+dN,:)-simVel(t-tau,:)))/sum(r.^2);
        lambdaDot=cross(r, (crickvel-handvel))/sum(r.^2);

        %Calculates commanded acceleration under delayed PN guidance law
        birdAcc2=N*cross(lambdaDot,simVel(t-1,:));     
        % birdAcc2 = N * cross( (lureVel(t-1,:)-simVel(t-1,:)), lambdaDot);  % 2D implementation
        
        %% Computes bird's velocity and position at time step t
    
        %Applies commanded acceleration to bird's simulated velocity
        simVel(t,:)=simVel(t-1,:) + (birdAcc+birdAcc2)/stepSize;
        
        %Adjusts simulated velocity of bird to match its measured speed
        simVel(t,:)=birdSpeed(t,:)*simVel(t,:)/norm(simVel(t,:));
        
        %Computes position of bird at time step t
        simPos(t,:)=simPos(t-1,:)+simVel(t,:)/stepSize;              
  end

return;


%%****** OLDER CODE DID NOT ACCOUNT FOR NEGATIVE TIMES 

% %Copyright. Graham K. Taylor & Caroline H. Brighton, University of Oxford, November 2018.
% 
% function [PPx,PPy,r2,N] = FindBestPNP(Hx,Hy,Cx,Cy,K,N,tau,lamb,minoff,taup)
%       
%       stepSize = 240;
%       tauMax = max(tau,minoff);
%       if (tauMax > (length(Hx)-1) )
%           PPx = Hx;
%           PPy = Hy;
%           r2 = NaN;
%           N = 0;
%           return;
%       end      
%       birdPos = [Hx,Hy,zeros(size(Hx))];
%       birdVel = [diff(Hx),diff(Hy),zeros(size(diff(Hx)))] * stepSize;
%       birdVel = [birdVel(1,:) ; birdVel];
%       birdSpeed = vecnorm(birdVel')';
%       lurePos = [Cx,Cy,zeros(size(Cx))];
%       lureVel = [diff(Cx),diff(Cy),zeros(size(diff(Cx)))] * stepSize;
%       lureVel = [lureVel(1,:) ; lureVel];
%       %****** call their routine
%       if isnan(tau)
%           input('what?');
%       end
%       simPos = zfitPNP(birdPos,birdVel,birdSpeed,lurePos,lureVel,stepSize,tauMax,K,N,tau,taup);
%       %*******
%       PPx = simPos(:,1);
%       PPy = simPos(:,2);
%       ex = tauMax:length(Hx);
%       Ptot = sum( (Hx(ex)-mean(Hx(ex))) .^ 2) + sum( (Hy(ex)-mean(Hy(ex))) .^ 2 );
%       Perr = sum( (Hx(ex) - PPx(ex)) .^ 2) + sum( (Hy(ex) - PPy(ex)) .^ 2);
%       Vtot = sum( diff(Hx(ex)) .^ 2) + sum( diff(Hy(ex)) .^ 2);
%       Verr = sum( (diff(Hx(ex)) - diff(PPx(ex))) .^ 2) + sum( (diff(Hy(ex)) - diff(PPy(ex))) .^ 2);    
%       %*******
%       N = (length(Hx)-tauMax);
%       r2 = ((1-lamb)*((Ptot-Perr)/Ptot) + lamb * ((Vtot-Verr)/Vtot));   % weight position and vel R2
%       %********
% return
% 
% 
% %Copyright. Graham K. Taylor & Caroline H. Brighton, University of Oxford, November 2018.
% function [simPos] = zfitPNP(birdPos,birdVel,birdSpeed,lurePos,lureVel,stepSize,tauMax,K,N,tau,taup)
% 
%   %Initializes matrices for bird's simulated position and velocity
%   simPos=birdPos*nan; simPos(1:tauMax,:)=birdPos(1:tauMax,:);
%   simVel=birdVel*nan; % simVel(1:tauMax,:)=birdVel(1:tauMax,:);
%   
%   dN = size(lurePos,1) - size(birdPos,1);  % cricket trace is appended (sub dN to sync times)
% 
%   %*** rest to straight line for history up to TauMax (reduce noise)
%   uvel = nanmean(birdVel(1:tauMax,:));
%   for t = 1:tauMax
%       simVel(t,:) = uvel;
%   end
%   %*********
%   
%   for t=tauMax+1:length(birdPos)
%     
%     %% Uses PP guidance law to compute bird's commanded acceleration at time step (t-1)
%       
%         %Calculates line-of-sight vector (r) lagged by tau
%         % r=lurePos(t-tau,:)-simPos(t-tau,:);
%         %**** Note, now r can be updated by future prediction "taup"
%         r = ((taup * lureVel(t-tau+dN)) + lurePos(t-tau+dN,:)) - simPos(t-tau,:);
%     
%         %Defines velocity vector (v) lagged by tau
%         v=simVel(t-tau,:);
%     
%         if (0)
%           %Calculates deviation angle (delta) lagged by tau
%           delta=acos(dot(r,v)/(norm(r)*norm(v)));
%     
%           %Creates unit vector normal to line of sight (r) and velocity (v)
%           n=cross(r,v)/(norm(r)*norm(v)); n=n/norm(n); 
%     
%           %Calculates commanded acceleration under delayed PP guidance law
%           birdAcc=-K*delta*cross(n,simVel(t-1,:));
%      
%         else
%           % use direct 2D version of pure pursuit (not angle based steering)
%           birdAcc = K * r * stepSize;
%         end
%         
%         %% Uses PN guidance law to compute bird's commanded acceleration at time (t-1)
%         %Calculates line-of-sight vector (r) lagged by tau
% 
%         %Calculates line-of-sight rate (lambdaDot) at time (t-delay)
%         lambdaDot=cross(r, (lureVel(t-tau+dN,:)-simVel(t-tau,:)))/sum(r.^2);
% 
%         %Calculates commanded acceleration under delayed PN guidance law
%         birdAcc2=N*cross(lambdaDot,simVel(t-1,:));
%         
%         %% Computes bird's velocity and position at time step t
%         %Applies commanded acceleration to bird's simulated velocity
%         simVel(t,:)=simVel(t-1,:) + (birdAcc+birdAcc2)/stepSize;
%         
%         %Adjusts simulated velocity of bird to match its measured speed
%         simVel(t,:)=birdSpeed(t,:)*simVel(t,:)/norm(simVel(t,:));
%         
%         %Computes position of bird at time step t
%         simPos(t,:)=simPos(t-1,:)+simVel(t,:)/stepSize;              
%   end
% 
% return;
% 
