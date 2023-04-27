
function pursuit_model_pred_delay_fpm()

% Generates Figure 3I, Figure 4F from
% Shaw,L, Wang KH, Mitchell, J (2024) Fast Prediction in Marmoset Reach-to-Grasp Movements for Dynamic Prey.
%
% This script compares pure pursuit, proportional navigation, and a mixed
% pursuit simulation performance at different visuomotor delays. Furthermore,
% the mixed pursuit simulation is compared to a predictive pursuit
% simulation.
% 
% Reaching data structure marmo_reach_model.mat available at
% https://doi.org/10.5281/zenodo.7869286
%
% Luke Shaw, Kuan Hong Wang, and Jude Mitchell 4/2023
% Matlab R2022b
%
% Pursuit simulation scripts modified from:
% Brighton, C. H. and G. K. Taylor (2019). "Hawks steer attacks using a guidance system tuned for close pursuit of erratically manoeuvring targets." Nat Commun 10(1): 2462.

%%

PreloadFile = 'ReachSimPredDelay.mat';

if ~exist(PreloadFile)
    
  AutoReg = 0;  % if 1, correct position each step (fit velocity error)
  MarkTau = 17;
  shouldVis = 0; %1=yea vis 0=no vis

  
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

  lamb = 1.0;  % 0 - weight position fit, 1 weight velocity fit, all 80ms delay 
                 % using velocity error, R2 is PP: 89.1(1.5)
                 %                             PN: 93.3(1.3)
                 %                             PP2: 93.6(1.1)
                 %                             PN2: 90.2(1.9)
  if (AutoReg == 1)
      lamb = 1.0;  % only makes sense to weight velocity errors
  end
  minoff = 18; % force identical following this number of steps
  
  %******* fit entire dataset with one parameterization
  %*********
  Gmin = 0.1;
  Gmax = 50;
  Gn = 40; % set higher for final run 50;
  GSAMP = exp( log(Gmin):((log(Gmax)-log(Gmin))/Gn):log(Gmax) );
  % TauPset = -0.5:0.05:0.5; 
  TauPset = -0.50:0.02:0.50;
  TauN = length(TauPset);
  % VMset = 1:4:40; 
  VMset = 1:2:40;
  VMN = length(VMset);
 
  %*********** everything you need to know for making plots, stored here
  Keep = struct;
  Keep.TauPset = TauPset;
  Keep.TauN = TauN;
  Keep.VMset = VMset;
  Keep.VMN = VMN;
  Keep.Gn = Gn;
  Keep.GSAMP = GSAMP;
  Keep.R2 = NaN(TauN,VMN,4);  % three models, PP, PN, PP(pred)
  Keep.SR2 = NaN(TauN,VMN,4);  % three models, PP, PN, PP(pred)
  Keep.Gains = NaN(TauN,VMN,4);  % gain parameters for each model
  Keep.Best = cell(TauN,VMN,4); % store trial by trial R2 and N
 %********* use capital to find the minimum for showing example (DTAU)
 BestPP_R2 = 0;
 BestPP_gain = NaN;
 BestPN_R2 = 0;
 BestPN_gain = NaN;
 BestPP2_R2 = 0;
 BestPP2_taup = 0;
 BestPP2_gain = NaN;
 BestPN2_R2 = 0;
 BestPN2_gain = NaN;
 BestPN2_gain2 = NaN;
 %*******
 for dtau = 1:VMN
  DTau = VMset(dtau);
  disp(['Stepping VM Tau ',num2str(DTau)]);   
  for tp = 1:TauN
    gp = (1+floor((tp-1)*Gn/TauN));  
    gain2 = GSAMP(gp);
    taup = TauPset(tp);  
    disp(['Stepping taup ',num2str(taup)]);  
    %****** for each of the three models, find best gain and error
    bestPP_R2 = 0;
    bestPP_gain = NaN;
    bestPPtrials = [];
    bestPN_R2 = 0;
    bestPN_gain = NaN;
    bestPNtrials = [];
    bestPP2_R2 = 0;
    bestPP2_gain = NaN;
    bestPP2trials = [];
    bestPN2_R2 = 0;
    bestPN2_gain = NaN;
    bestPN2_gain2 = NaN;
    bestPN2trials = [];
    %****************
    for gain = GSAMP % Brighton + Taylor, K or N term
         % disp(['Stepping ',num2str(gain)]);  
         %********************        
         PPsum = 0; sumPPn = 0;
         PNsum = 0; sumPNn = 0;
         PP2sum = 0; sumPP2n = 0;
         PN2sum = 0; sumPN2n = 0;
         BestTrialSet = cell(1,4);
         for i = 1:nReach % number of reaches
           Hx= smoothdata(model.x.hand{i}(1:end,:),'gaussian',smw);
           Hy= smoothdata(model.y.hand{i}(1:end,:),'gaussian',smw);
           zCx = vertcat(model.x.cricketexfull{i},model.x.cricket{i});
           zCy = vertcat(model.y.cricketexfull{i},model.y.cricket{i});
           zCx = smoothdata(zCx,'gaussian',smw);
           zCy = smoothdata(zCy,'gaussian',smw);
           Cx= smoothdata(model.x.cricket{i}(1:end,:),'gaussian',smw);
           Cy= smoothdata(model.y.cricket{i}(1:end,:),'gaussian',smw);
           
           %******allow trial by trial start direction search (to account
           %****** for arbitrary start posture with each trial)
           mnoff = minoff;
           % mnoff = DTau;
           if (tp == 1)
              [PPx,PPy,PPr2,PPn] = FindBestPNP(Hx,Hy,zCx,zCy,gain,0,DTau,lamb,mnoff,0,0,AutoReg);
              [PNx,PNy,PNr2,PNn] = FindBestPNP(Hx,Hy,zCx,zCy,0,gain,DTau,lamb,mnoff,0,0,AutoReg);       
              %****
              if ~isnan(PPr2)
                PPsum = PPsum + (PPn * PPr2);
                sumPPn = sumPPn + PPn;
              end
              if ~isnan(PNr2)
                PNsum = PNsum + (PNn * PNr2);
                sumPNn = sumPNn + PNn;
              end
              BestTrialSet{1} = [BestTrialSet{1} ; [PPr2, PPn]];
              BestTrialSet{2} = [BestTrialSet{2} ; [PNr2, PNn]];
           end
           %******* prediction variances
           [PP2x,PP2y,PP2r2,PP2n] = FindBestPNP(Hx,Hy,zCx,zCy,gain,0,DTau,lamb,mnoff,taup,1,AutoReg);  % 2D prediction model
           [PN2x,PN2y,PN2r2,PN2n] = FindBestPNP(Hx,Hy,zCx,zCy,gain,gain2,DTau,lamb,mnoff,0,0,AutoReg);
           %****
           BestTrialSet{3} = [BestTrialSet{3} ; [PP2r2, PP2n]];
           if ~isnan(PP2r2)
                PP2sum = PP2sum + (PP2n * PP2r2);
                sumPP2n = sumPP2n + PP2n;
           end
           BestTrialSet{4} = [BestTrialSet{4} ; [PN2r2, PN2n]];
           if ~isnan(PN2r2)
                PN2sum = PN2sum + (PN2n * PN2r2);
                sumPN2n = sumPN2n + PN2n;
           end
           %*******
         end % nReach, best fit over total set
         
         if (tp == 1)
           PPsum = (PPsum/sumPPn);  % norm over trials, total R2
           PNsum = (PNsum/sumPNn);  % norm over trials, total R2
           %*********************
           if (DTau == MarkTau)
              if (PPsum > BestPP_R2)
                 BestPP_gain = gain;
                 BestPP_R2 = PPsum;
              end
           end
           if (PPsum > bestPP_R2)
              bestPP_gain = gain;
              bestPP_R2 = PPsum;
              bestPPtrials = BestTrialSet{1};
           end
           %*****************
           if (DTau == MarkTau)
              if (PNsum > BestPN_R2)
                 BestPN_gain = gain;
                 BestPN_R2 = PNsum;
              end
           end
           if (PNsum > bestPN_R2)
              bestPN_gain = gain;
              bestPN_R2 = PNsum;
              bestPNtrials = BestTrialSet{2};
           end
         end
         %********** compute for prediction variants
         PP2sum = (PP2sum/sumPP2n);  % norm over trials, total R2
         PN2sum = (PN2sum/sumPN2n);  % norm over trials, total R2
         %*********************
         if (DTau == MarkTau)
            if (PP2sum > BestPP2_R2)
                BestPP2_R2 = PP2sum;
                BestPP2_gain = gain;
                BestPP2_taup = taup;
            end
         end
         if (PP2sum > bestPP2_R2)
              bestPP2_gain = gain;
              bestPP2_taup = taup;
              bestPP2_R2 = PP2sum;
              bestPP2trials = BestTrialSet{3};
         end
         %*********************
         if (DTau == MarkTau)
            if (PN2sum > BestPN2_R2)
                BestPN2_R2 = PN2sum;
                BestPN2_gain = gain;
                BestPN2_gain2 = gain2;
            end
         end
         if (PN2sum > bestPN2_R2)
              bestPN2_gain = gain;
              bestPN2_gain2 = gain2;
              bestPN2_R2 = PN2sum;
              bestPN2trials = BestTrialSet{4};
         end
         %*******
    end % Gain, fit per DTau and TauP
    %***** first set on the first two models
    if (tp == 1)
       Keep.R2(tp,dtau,1) = bestPP_R2; 
       Keep.Gain(tp,dtau,1) = bestPP_gain;
       Keep.Best{tp,dtau,1} = bestPPtrials;
       if ~isempty(bestPPtrials)
          N2 = sum(~isnan(bestPPtrials(:,1)));
          if (N2)
            Keep.SR2(tp,dtau,1) = nanstd(bestPPtrials(:,1))/sqrt(N2);
          end
       end
       %*********
       Keep.R2(tp,dtau,2) = bestPN_R2; 
       Keep.Gain(tp,dtau,2) = bestPN_gain;
       Keep.Best{tp,dtau,2} = bestPNtrials;
       if ~isempty(bestPNtrials)
         N2 = sum(~isnan(bestPNtrials(:,1)));
         if (N2)
            Keep.SR2(tp,dtau,2) = nanstd(bestPNtrials(:,1))/sqrt(N2);
         end
       end
       %********* no need to assign others
    else
       Keep.R2(tp,dtau,1) = Keep.R2(1,dtau,1);
       Keep.R2(tp,dtau,2) = Keep.R2(1,dtau,2);
       Keep.SR2(tp,dtau,1) = Keep.SR2(1,dtau,1);
       Keep.SR2(tp,dtau,2) = Keep.SR2(1,dtau,2);
       Keep.Gain(tp,dtau,1) = Keep.Gain(1,dtau,1);
       Keep.Gain(tp,dtau,2) = Keep.Gain(1,dtau,2);
       Keep.Best{tp,dtau,1} = Keep.Best{1,dtau,1};
       Keep.Best{tp,dtau,2} = Keep.Best{1,dtau,2};
    end
    Keep.R2(tp,dtau,3) = bestPP2_R2; 
    Keep.Gain(tp,dtau,3) = bestPP2_gain;
    Keep.Best{tp,dtau,3} = bestPP2trials;
    if ~isempty(bestPP2trials)
        N2 = sum(~isnan(bestPP2trials(:,1)));
        if (N2)
            Keep.SR2(tp,dtau,3) = nanstd(bestPP2trials(:,1))/sqrt(N2);
        end
    end
    %************
    Keep.R2(tp,dtau,4) = bestPN2_R2; 
    Keep.Gain(tp,dtau,4) = bestPN2_gain;
    Keep.Best{tp,dtau,4} = bestPN2trials;
    if ~isempty(bestPN2trials)
        N2 = sum(~isnan(bestPN2trials(:,1)));
        if (N2)
            Keep.SR2(tp,dtau,4) = nanstd(bestPN2trials(:,1))/sqrt(N2);
        end
    end
    %*************
  end  % TauP
 end % DTau (VM delay in this sim)

 if (shouldVis == 1)
    for i = 1:nReach % number of reaches
      %**********
      Hx= smoothdata(model.x.hand{i}(1:end,:),'gaussian',smw);
      Hy= smoothdata(model.y.hand{i}(1:end,:),'gaussian',smw);
      zCx = vertcat(model.x.cricketexfull{i},model.x.cricket{i});
      zCy = vertcat(model.y.cricketexfull{i},model.y.cricket{i});
      zCx = smoothdata(zCx,'gaussian',smw);
      zCy = smoothdata(zCy,'gaussian',smw);
      Cx= smoothdata(model.x.cricket{i}(1:end,:),'gaussian',smw);
      Cy= smoothdata(model.y.cricket{i}(1:end,:),'gaussian',smw);
      %******
      DTau = MarkTau;
      %*********
      [PPx,PPy,PPr2,PPn] = FindBestPNP(Hx,Hy,Cx,Cy,BestPP_gain,0,DTau,lamb,mnoff,0,0,AutoReg);
      [PNx,PNy,PNr2,PNn] = FindBestPNP(Hx,Hy,Cx,Cy,0,BestPN_gain,DTau,lamb,mnoff,0,0,AutoReg);       
      %******* prediction variances
      [PP2x,PP2y,PP2r2,PP2n] = FindBestPNP(Hx,Hy,Cx,Cy,BestPP2_gain,0,DTau,lamb,mnoff,BestPP2_taup,1,AutoReg);  % 2D prediction model
      [PN2x,PN2y,PN2r2,PN2n] = FindBestPNP(Hx,Hy,Cx,Cy,BestPN2_gain,BestPN2_gain2,DTau,lamb,mnoff,0,0,AutoReg);
      %*************
      if (1)
         hf = figure(2); hold off;
         set(hf,'Position',[100 200 900 500]);
         skip = 4;
         rcolo = [1,0.5,0.5];
         bcolo = [0,0.5,0.5];
         gcolo = [0.7,0.7,0.7];
         for tt = (1+skip):skip:size(PP2x,1)
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
              plot([PP2x(tt),PP2x(tt-skip)],[PP2y(tt),PP2y(tt-skip)],'r-');
              %**********
              plot([PN2x(tt),PN2x(tt-skip)],[PN2y(tt),PN2y(tt-skip)],'m-');
              %*******
         end
         axis([0 20 -2 10]);
         ht = title(sprintf(['T(%d) PP(%4.3f,%3.1f,%d) PN(%4.3f,%3.1f,%d)',...
                           ' PF(%4.3f,%3.1f,%d) PND(%4.3f,%3.1f,%3.1f,%d)'],...
                      i,PPr2,BestPP_gain,DTau,...
                        PNr2,BestPN_gain,DTau,...
                        PP2r2,BestPP2_gain,DTau,...
                        PN2r2,BestPN2_gain,BestPN2_gain2,DTau));
         set(ht,'Fontsize',10);
         [PPr2,PNr2,PP2r2,PN2r2]
         input('check');
      end
    end
 end
 %***************
 disp(['Saving preloaded data to ',PreloadFile]);
 save(PreloadFile,'Keep');
 disp('Finished saving');
else
  disp(['Loading preloaded data from ',PreloadFile]);
  load(PreloadFile);  
  disp('Loaded');
end


%******** resort the mixed PNP model, find best gain2 and copy
VMN = size(Keep.R2,2);
TauN = size(Keep.R2,1);
for k = 1:VMN
    sq2 = squeeze(Keep.R2(:,k,4));
    maxo = max(sq2);
    zz = find( sq2 == maxo);
    tp = zz(1);
    for kk = 1:TauN
       if (kk ~= tp)
          Keep.R2(kk,k,4) = Keep.R2(tp,k,4);
          Keep.Gain(kk,k,4) = Keep.Gain(tp,k,4);
          Keep.Best{kk,k,4} = Keep.Best{tp,k,4};
          Keep.SR2(kk,k,4) = Keep.SR2(tp,k,4);
       end
    end
end

%**** preliminary, make some fast plots *********
DTAU = 80;
tscal = (1000/240);
if (0)
    hf = figure;
    set(hf,'Position',[100 100 1800 500]);
    zmaxo = max(max(max(Keep.R2)));
    % zmino = min(min(min(Keep.R2)));
    zmino = 0.88;
    dmino = 0.00;
    for k = 1:6
      subplot('Position',[(0.05+((k-1)*0.16)) 0.15 0.125 0.70]);
      if (k < 4)
         imo = squeeze(Keep.R2(:,:,k+1));
         imagesc(Keep.VMset*tscal,Keep.TauPset,imo,[zmino zmaxo]); hold on;
      else
         if (k == 4)
             imo = squeeze(Keep.R2(:,:,3)-Keep.R2(:,:,2));
             imagesc(Keep.VMset*tscal,Keep.TauPset,imo,[dmino 0.07]); hold on;    
         else
            if (k == 5)
              imo = squeeze(Keep.R2(:,:,3)-Keep.R2(:,:,4));
              imagesc(Keep.VMset*tscal,Keep.TauPset,imo,[dmino max(max(imo))]); hold on;
            else
              imo = squeeze(Keep.R2(:,:,3)-Keep.R2(:,:,4));
              simo = squeeze( 0.5*(Keep.SR2(:,:,3)+Keep.SR2(:,:,4)));
              zimo = imo ./ simo;
              imagesc(Keep.VMset*tscal,Keep.TauPset,zimo,[0 5]); hold on;       
            end
         end
      end
      %***** find the peak at 80 ms delay
      vms = Keep.VMset*tscal;
      disto = (Keep.VMset*tscal - DTAU).^2;
      zm = find( disto == min(disto) );
      vm = zm(1);
      maxo = max(max(imo(:,vm)));
      zz = find( imo == maxo);
      a = 1+floor((zz(1)-1)/Keep.TauN);
      b = mod((zz(1)-1),Keep.TauN)+1;
      if (k >= 4)
          a5 = a;
          b5 = b;
          disp('Peak fit taup and VM delay');
          [Keep.TauPset(b),Keep.VMset(a)*tscal]
      end
      if (k == 5)  % do stats compare
          bestPP2 = Keep.Best{b5,a5,3};  % errors per trial at peak
          bestPN = Keep.Best{b5,a5,2};
          zz = find(~isnan(bestPP2(:,1)));
          % p1 = signrank(bestPP2(zz,1),bestPN(zz,1));  % PP2 vs PN model
          [h1,p1] = ttest(bestPP2(zz,1),bestPN(zz,1));
          disp('Sig dif PP2 vs PN?');
          [(sum( bestPP2(zz,1) .* bestPP2(zz,2))/sum(bestPP2(zz,2))),...
           (sum( bestPN(zz,1) .* bestPN(zz,2))/sum(bestPN(zz,2)))]
          p1

          %*******
          bestPP2 = Keep.Best{b5,a5,3};  % errors per trial at peak
          bestPN2 = Keep.Best{b5,a5,4};
          [b5,a5]
          if ~isempty(bestPN2)
            zz = find(~isnan(bestPN2(:,1)));
            % p2 = signrank(bestPP2(zz,1),bestPN2(zz,1));  % PP2 vs PN model
            [h2,p2] = ttest(bestPP2(zz,1),bestPN2(zz,1));
            disp('Sig dif PP2 vs PN2?');
            [(sum( bestPP2(zz,1) .* bestPP2(zz,2))/sum(bestPP2(zz,2))),...
             (sum( bestPN2(zz,1) .* bestPN2(zz,2))/sum(bestPN2(zz,2)))]
            p2
          end

      end
      plot(Keep.VMset(a)*tscal,Keep.TauPset(b),'k+');
      %*******
      colorbar;
      xlabel('VM delay');
      ylabel('taup ');
      switch k
          case 1  
              title('PN err');
          case 2
              title('Pred err');
          case 3
              title('Mix err');
          case 4
              title('Pred - PN');
          case 5 
              title('Pred - Mix');
          case 6
              title('Zscore (Pred-Mix)');
       end
    end
end
%**** look at difference maps


%*** make final sub-panel to Figure 3, break-down by visuo-motor delay
if (1)
    includePNP = 1;
    includePF = 1;
    PPcolo = [0,0,0];
    PNcolo = [0,0,1];
    PFcolo = [0.6,0.4,0];
    PN2colo = [0.8,0.0,0.4];
    %****
    hf = figure;
    set(hf,'Position',[100 100 600 500]);
    xtime = Keep.VMset * tscal;
    uPP = squeeze(Keep.R2(1,:,1));
    sPP = squeeze(Keep.SR2(1,:,1));
    uPN = squeeze(Keep.R2(1,:,2));
    sPN = squeeze(Keep.SR2(1,:,2));
    uPF = [];
    sPF = [];  % must search for max over each taup
    uPN2 = [];
    sPN2 = [];
    for k = 1:size(Keep.R2,2)
       slice = squeeze(Keep.R2(:,k,3));
       maxo = max(slice);
       zz = find( slice == maxo);
       uPF = [uPF  Keep.R2(zz(1),k,3)];
       sPF = [sPF  Keep.SR2(zz(1),k,3)];
       %**********
       slice = squeeze(Keep.R2(:,k,4));
       maxo = max(slice);
       zz = find( slice == maxo);
       uPN2 = [uPN2  Keep.R2(zz(1),k,4)];
       sPN2 = [sPN2  Keep.SR2(zz(1),k,4)];
    end
    %**** plot curves
    aa = [xtime fliplr(xtime)];
    bb = [(uPP+sPP) fliplr(uPP-sPP)];
    fill(aa,bb,PPcolo,'FaceAlpha',0.3,'Linestyle','none'); hold on;
    plot(xtime,uPP,'k-','LineWidth',2,'Color',PPcolo);
    bb = [(uPN+sPN) fliplr(uPN-sPN)];
    fill(aa,bb,PNcolo,'FaceAlpha',0.3,'Linestyle','none'); hold on;
    plot(xtime,uPN,'k-','LineWidth',2,'Color',PNcolo);
    pnbest = max(uPN);
    plot(xtime,(pnbest*ones(size(xtime))),'b:','Color',PNcolo,'Linewidth',2);
    %***** leave out PD controller for now
    if (includePF == 1) 
      bb = [(uPF+sPF) fliplr(uPF-sPF)];
      fill(aa,bb,PFcolo,'FaceAlpha',0.3,'Linestyle','none'); hold on;
      plot(xtime,uPF,'k-','LineWidth',2,'Color',PFcolo);  
      text(80,0.960,'Predictive Pursuit','Fontsize',14,'Color',PFcolo);
    end
    if (includePNP == 1)
      bb = [(uPN2+sPN2) fliplr(uPN2-sPN2)];
      fill(aa,bb,PN2colo,'FaceAlpha',0.3,'Linestyle','none'); hold on;
      plot(xtime,uPN2,'k-','LineWidth',2,'Color',PN2colo);
      text(80,0.975,'Mixed PP + PN','Fontsize',14,'Color',PN2colo);
    end
    plot([DTAU,DTAU],[0.78,0.95],'k:');  % show 80 ms delay
    xlabel('Visuomotor Delay (ms)');
    ylabel('Goodness of Fit');
    set(gca,'Fontsize',14);
    axis([4 160 0.78 1.02]);
 end

%*** make final sub-panel to Figure 3, break-down by visuo-motor delay
if (1)
    includePNP = 1;
    includePF = 1;
    PFcolo = [0.6,0.4,0];
    PN2colo = [0.8,0.0,0.4];
    %****
    hf = figure;
    set(hf,'Position',[100 100 600 300]);
    subplot('Position',[0.15 0.20 0.70 0.70]);
    xtime = Keep.VMset * tscal;
    uPF = [];
    sPF = [];  % must search for max over each taup
    uPN2 = [];
    sPN2 = [];
    for k = 1:size(Keep.R2,2)
       slice = squeeze(Keep.R2(:,k,3));
       maxo = max(slice);
       zz = find( slice == maxo);
       uPF = [uPF  Keep.R2(zz(1),k,3)];
       sPF = [sPF  Keep.SR2(zz(1),k,3)];
       %**********
       slice = squeeze(Keep.R2(:,k,4));
       maxo = max(slice);
       zz = find( slice == maxo);
       uPN2 = [uPN2  Keep.R2(zz(1),k,4)];
       sPN2 = [sPN2  Keep.SR2(zz(1),k,4)];
    end
    %**** plot curves
    colororder([[0.6,0.0,0.3];[0,0,0]]);
    yyaxis left;
    aa = [xtime fliplr(xtime)];
    if (includePF == 1) 
      bb = [(uPF+sPF) fliplr(uPF-sPF)];
      fill(aa,bb,PFcolo,'FaceAlpha',0.3,'Linestyle','none'); hold on;
      plot(xtime,uPF,'k-','LineWidth',2,'Color',PFcolo);  
      text(10,0.955,'Predictive Pursuit','Fontsize',14,'Color',PFcolo);
    end
    if (includePNP == 1)
      bb = [(uPN2+sPN2) fliplr(uPN2-sPN2)];
      fill(aa,bb,PN2colo,'FaceAlpha',0.3,'Linestyle','none'); hold on;
      plot(xtime,uPN2,'k-','LineWidth',2,'Color',PN2colo);
      text(10,0.970,'Mixed Pursuit','Fontsize',14,'Color',PN2colo);
    end
    % plot([DTAU,DTAU],[0.85,0.98],'k--','Linewidth',1.5);  % show 80 ms delay
    ymin = 0.85;
    ymax = 0.98;
    %gaa = [(DTAU-5) (DTAU-5) (DTAU+5) (DTAU+5) (DTAU-5)];
    %gbb = [ymin ymax ymax ymin ymin];
    % fill(gaa,gbb,[0,0,0],'FaceAlpha',0.2,'Linestyle','none');
    xlabel('Visuomotor Delay (ms)');
    ylabel('Goodness of Fit');
    set(gca,'Fontsize',14);
    axis([4 160 0.85 0.98]);
    xticks([0 40 80 120]);
    yticks([0.85 0.90 0.95]);
    %***********
    yyaxis right;
    uu = uPF - uPN2;
    su = 0.5*(sPF+sPN2);
    bb = [(uu+su) fliplr(uu-su)];
    fill(aa,bb,[0,0,0],'FaceAlpha',0.1,'Linestyle','none');
    plot(xtime,uu,'k-','Linewidth',2); hold on;
    plot(xtime,uu+su,'k-','Linewidth',1);
    plot(xtime,uu-su,'k-','Linewidth',1);
    plot(xtime,zeros(size(xtime)),'k-','Linewidth',2);
    ylabel('Predictive - Mixed');
    yticks([0.00 0.04]); % yticks([0 0.02 0.04]);
    set(gca,'Fontsize',14);
    set(gca,'Linewidth',1.5);
    axis([4 160 -0.005 0.045]);
    disto = (xtime-DTAU).^2;
    dd = find( disto == min(disto) );
    div = 0.001;
    plot([DTAU,DTAU],[div,((uPF(dd)-uPN2(dd))-div)],'k-','Linewidth',3);
    text((DTAU+5),(0.35*(uPF(dd)-uPN2(dd))),'*','Fontsize',24);
    %********
    if (1)  % do stats test on model diff at 80ms
      imo = Keep.R2(:,:,3);  % PF pred model
      maxo = max(max(imo(:,dd)));
      zz = find( imo == maxo);
      a = 1+floor((zz(1)-1)/Keep.TauN);
      b = mod((zz(1)-1),Keep.TauN)+1;
      disp('Peak fit taup and VM delay');
      [Keep.TauPset(b),Keep.VMset(a)*tscal]
      %****
      bestPP2 = Keep.Best{b,a,3};  % errors per trial at peak
      bestPN2 = Keep.Best{b,a,4};
      zz = find(~isnan(bestPP2(:,1)));
      psign = signrank(bestPP2(zz,1),bestPN2(zz,1));  % PP2 vs PN model
      [h1,pttest] = ttest(bestPP2(zz,1),bestPN2(zz,1));
      disp(sprintf('Sig dif at %d: Predictive - Mixed',Keep.VMset(a)*tscal));
      [(sum( bestPP2(zz,1) .* bestPP2(zz,2))/sum(bestPP2(zz,2))),...
       (sum( bestPN2(zz,1) .* bestPN2(zz,2))/sum(bestPN2(zz,2)))]
      psign
      pttest
      title(sprintf('Signtest(p=%6.4f)  t-test(p=%6.4f)',psign,pttest),'Fontsize',12);
    end
 end


return;


%********* BELOW, TRYING TO INCORPORATE USING EXTENDED CRICK VELOCITIES 
%Copyright. Graham K. Taylor & Caroline H. Brighton, University of Oxford, November 2018.
function [PPx,PPy,r2,N] = FindBestPNP(Hx,Hy,Cx,Cy,K,N,tau,lamb,minoff,taup,use2D,AutoReg)
      
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
      simPos = zfitPNP(birdPos,birdVel,birdSpeed,lurePos,lureVel,stepSize,tauMax,K,N,tau,taup,use2D,AutoReg);
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
function [simPos] = zfitPNP(birdPos,birdVel,birdSpeed,lurePos,lureVel,stepSize,tauMax,K,N,tau,taup,use2D,AutoReg)

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
        %**** Note, now r can be updated by future prediction "taup"
        %**** and can implement reverse approximation if cricket position missing
        if ((t-tau+dN) < 1)
            db = t-tau+dN-1;
            crickpos = lurePos(1,:) + ((db/stepSize) * cvel(1,:));
            crickvel = cvel(1,:);
        else
            crickpos = lurePos(t-tau+dN,:);
            crickvel = lureVel(t-tau+dN,:);
        end
        
        if (0)  % assume zero delay on hand position information (worse performance)
           handpos = simPos(t-1,:); 
           handvel = simVel(t-1,:);
        else 
           if ( (t-tau)<1)  % implement reverse extrapolate if negative time
              db = t-tau-1;
              handpos = simPos(1,:) + ((db/stepSize) * simVel(1,:));
              handvel = simVel(1,:);
           else
              if (AutoReg == 1)
                 handpos = birdPos(t-tau,:);  % use real position
                 handvel = birdVel(t-tau,:);    
              else
                 handpos = simPos(t-tau,:);
                 handvel = simVel(t-tau,:);
              end
           end
        end
        
        %******* prediction from cricket position
        lurePred = (taup * crickvel) + crickpos;
        r = lurePred - handpos;  %effectively using hand position a delay
         
        %Defines velocity vector (v) lagged by tau
        v=handvel;  % simVel(t-tau,:);
    
        if ~use2D  % Steering pure pursuit does much worse in this system
          %Calculates deviation angle (delta) lagged by tau
          delta=acos(dot(r,v)/(norm(r)*norm(v)));
    
          %Creates unit vector normal to line of sight (r) and velocity (v)
          n=cross(r,v)/(norm(r)*norm(v)); n=n/norm(n); 
    
          %Calculates commanded acceleration under delayed PP guidance law
          birdAcc=-K*delta*cross(n,simVel(t-1,:));
      
        else  % hold out 2D pure pursuit ... that is a non-steering strategy
            
          %**** translate into 2D pursuit vector (forget steering)
          birdAcc = K * r * stepSize;
              
        end
        
        % yyaxis left or right
        
        %% Uses PN guidance law to compute bird's commanded acceleration at time (t-1)
        %Calculates line-of-sight vector (r) lagged by tau
        %Calculates line-of-sight rate (lambdaDot) at time (t-delay)
        % lambdaDot=cross(r, (lureVel(t-tau+dN,:)-simVel(t-tau,:)))/sum(r.^2);
        lambdaDot=cross(r, (crickvel-handvel))/sum(r.^2);

        %Calculates commanded acceleration under delayed PN guidance law
        birdAcc2=N*cross(lambdaDot,simVel(t-1,:));     
        % birdAcc2 = N * cross( (lureVel(t-1,:)-simVel(t-1,:)), lambdaDot);
        
        %% Computes bird's velocity and position at time step t
    
        %Applies commanded acceleration to bird's simulated velocity
        simVel(t,:)=simVel(t-1,:) + (birdAcc+birdAcc2)/stepSize;
        
        %Adjusts simulated velocity of bird to match its measured speed
        simVel(t,:)=birdSpeed(t,:)*simVel(t,:)/norm(simVel(t,:));
        
        %Computes position of bird at time step t
        %if (AutoReg == 1)
        %    simPos(t,:) = birdPos(t-1)+simVel(t,:)/stepSize;
        %else
            simPos(t,:) = simPos(t-1,:)+simVel(t,:)/stepSize;
        %end
  end

return;

