function estimate_delay_fpm(model,nSEM,nJack)

% Generates Figure 4A, 4B from  
% Shaw,L, Wang KH, Mitchell, J (2024) Fast Prediction in Marmoset Reach-to-Grasp Movements for Dynamic Prey.

% inputs:  
%   model   Reaching data structure marmo_reach_model.mat available at
%               https://doi.org/10.5281/zenodo.7869286
%   nSEM    SEM multiplier for plotting
%   nJack   number of jack knifes to perform on STA if 0 then confidence
%               bounds of regress function implemented
%
% This function estimates visumotor delay in marmoset reaches using two methods. 
% The first uses a filter that finds jumps in cricket velocity and calculates 
% lateral hand velocity locked to those events. The second is a simple
% linear regression model with a delay term. See Star Methods of above
% research publication.
% Jude Mitchell and Luke Shaw 4/2023
%
% Reaching data structure marmo_reach_model.mat available at
% https://doi.org/10.5281/zenodo.7869286

%% 

LT = length(model.x.hand);
OFPS = 240;  % frames per sec of video
FPS = 240;  % resample data to reduce auto-correlations
Lagtime = 25;  % lag for auto-regressive filter
MinT = 400; % in ms, min +/- lags
SGO = 3;
stype = 'sgolay';
%**** Use -25 to 150 as range in all plots
BigLag = floor((150+6.25)*FPS/1000);  % must be even, at 40hz this is 200 ms
AntLag = -floor((25+6.25)*FPS/1000); % provide some negative times (check correlation)
BMD = 3; % subsample BigX for STA before computation 
%*******
FP = 0.05;   % note, any earlier frames of hand would help this analysis
FPWIN = floor( FPS * 0.10);  % only flag events within +/- this period
FPWINA = floor( FPS * 0.075);
FPWINB = floor( FPS * 0.15);
JACKN = nJack; % if >1 , perform JACKNIFE ON STA
SUBTRIALMEANS = 0;  % doesnt change results (STA not due to trial-by-trial corr

%********* Acceleration filter (apply to abs velocity trace)
ACCSIG = floor(20*FPS/1000);  % gauss sig is 20ms, +/- 2 sig
AccFilt = zeros(1,(ACCSIG*4)+1);
SAccFilt = zeros(1,(ACCSIG*4)+1);
for ak = 1:((ACCSIG*4)+1)
    val = (ak-((ACCSIG*2)+1));
    if (val >= 0)
      AccFilt(1,ak) = exp(-0.5*(val/ACCSIG)^2);
      SAccFilt(1,ak) = AccFilt(1,ak);
    else
      AccFilt(1,ak) = -exp(-0.5*(val/ACCSIG)^2);
    end
end
AccFilt = AccFilt/norm(AccFilt);
%***********
withautoreg = 0;  % in linear regression include h(t-1)
PlotIt = 0;
%***** Build a reverse correlation model with time lags
BigX = [];
BigY = [];

%***** these trials do not have much cricket motion ... seem likely
%***** to just add noise to the analyses
MinC2 = floor((MinT*FPS/1000)/2);
CrossTime = (-MinC2:MinC2) * (1000/FPS);
MN = 1+(2*MinC2);
LAGC = 1+floor((Lagtime*FPS/1000));  % number of coefficients to invert autocorr
HandCross = [];
CrickCross = [];
HandCrick = [];
%***
HandCross3 = [];
CrickCross3 = [];
HandCrick3 = [];
%******
HandYule = [];
CrickYule = [];
%**** grand means
TOTclatvel = [];
TOThlatvel = [];
TOTclatacc = [];
TOThlatacc = [];
%% ******** computations to implement for lateral velocity
%******** and augment model information (will plot them below)
for tr = 1:LT
    %********* concatenate front part of criket
    if (1)
        cxpos = [model.x.cricketexfull{tr} ; model.x.cricket{tr}];
        cypos = [model.y.cricketexfull{tr} ; model.y.cricket{tr}];
        CN = length(model.x.cricketexfull{tr});
        NN = length(model.x.cricket{tr});
        hxpos = nan(size(cxpos));
        hypos = nan(size(cypos));
        hxpos((CN+1):(CN+NN)) = model.x.hand{tr};
        hypos((CN+1):(CN+NN)) = model.y.hand{tr};
    else
       cxpos = model.x.cricket{tr};
       cypos = model.y.cricket{tr};
       hxpos = model.x.hand{tr};
       hypos = model.y.hand{tr};
    end
    %*************
    model.velx.cricket{tr} = diff( cxpos ) * FPS;  % frames per sec
    model.velx.hand{tr} = diff( hxpos ) * FPS;  % frames per sec
    model.vely.cricket{tr} = diff( cypos ) * FPS;  % frames per sec
    model.vely.hand{tr} = diff( hypos ) * FPS;  % frames per sec
    %*********
    if (SGO) % smooth velocity
      %******* filter explicitly for increase in velocity
      movelxc = abs(model.velx.cricket{tr});
      movelyc = abs(model.vely.cricket{tr});
      movelxh = abs(model.velx.hand{tr});
      movelyh = abs(model.vely.hand{tr});
      %****
      Nac = length(AccFilt);
      Nac2 = floor(Nac/2)+1;
      %******* This filter looks for increases in abs speed
      Nco = length(movelxc);
      aco = conv(fliplr(AccFilt),movelxc');
      aso = conv(fliplr(SAccFilt),model.velx.cricket{tr}');
      xx = Nac2:(Nco+Nac2-1);
      model.accx.cricket{tr} = (sign( aso(xx) ) .* aco(xx))' * FPS;
      %*****
      Nco = length(movelyc);
      aco = conv(fliplr(AccFilt),movelyc');
      aso = conv(fliplr(SAccFilt),model.vely.cricket{tr}');
      xx = Nac2:(Nco+Nac2-1);
      model.accy.cricket{tr} = (sign( aso(xx) ) .* aco(xx))' * FPS;
      %********
      Nco = length(movelxh);
      aco = conv(fliplr(AccFilt),movelxh');
      aso = conv(fliplr(SAccFilt),model.velx.hand{tr}');
      xx = Nac2:(Nco+Nac2-1);
      model.accx.hand{tr} = (sign( aso(xx) ) .* aco(xx))' * FPS;
      %*************
      Nco = length(movelyh);
      aco = conv(fliplr(AccFilt),movelyh');
      aso = conv(fliplr(SAccFilt),model.vely.hand{tr}');
      xx = Nac2:(Nco+Nac2-1);
      model.accy.hand{tr} = (sign( aso(xx) ) .* aco(xx))' * FPS;
      
      %****
    end
    %**** update lateral velocity by projection position per time
    vecp = [(model.x.cricket{tr}(1)-model.x.hand{tr}(1)),...
        (model.y.cricket{tr}(1)-model.y.hand{tr}(1))];
    vecp = vecp/norm(vecp);
    latp = [-vecp(2) vecp(1)];  % orthogonal axis from reach
    clatvel = zeros(size(model.velx.cricket{tr}));
    hlatvel = zeros(size(model.velx.hand{tr}));
    if (SGO)
       clatacc = zeros(size(model.accx.cricket{tr}'));
       hlatacc = zeros(size(model.accx.hand{tr}'));
    end
    latp2=latp;
    for tt = 1:length(model.velx.hand{tr})
       %**** compute lateral position *****
       clatvel(tt) = (latp2(1) * model.velx.cricket{tr}(tt) + ...
                      latp2(2) * model.vely.cricket{tr}(tt));
       hlatvel(tt) = (latp2(1) * model.velx.hand{tr}(tt) + ...
                      latp2(2) * model.vely.hand{tr}(tt));
       %****** compute the ACC **************
       if (SGO) && (tt <= length(clatacc))
           clatacc(tt) = (latp2(1) * model.accx.cricket{tr}(tt) + ...
                          latp2(2) * model.accy.cricket{tr}(tt));
           hlatacc(tt) = (latp2(1) * model.accx.hand{tr}(tt) + ...
                          latp2(2) * model.accy.hand{tr}(tt));   
       end   
    end
    %******* eliminate trial by trial correlations (all zero mean)
    if (SUBTRIALMEANS)
       clatvel = clatvel - nanmean(clatvel);
       hlatvel = hlatvel - nanmean(hlatvel);
    end
    %**********
    if (OFPS == FPS)
       model.latvel.cricket{tr} = clatvel;
       model.latvel.hand{tr} = hlatvel; 
       clatvel2 = clatvel';
       hlatvel2 = hlatvel';
       NN = length(clatvel2);
       %********
       if (SGO)
          model.latacc.cricket{tr} = clatacc;
          model.latacc.hand{tr} = hlatacc; 
          aNN = length(hlatacc);
          %*******
          TOTclatacc = [TOTclatacc ; [clatacc' (tr*ones(size(clatacc'))) (1:aNN)']];
          TOThlatacc = [TOThlatacc ; [hlatacc' (tr*ones(size(hlatacc'))) (1:aNN)']];
          %*******          
       end
    else
       %***** down-sample to reduce correlations
       Sampo = floor( OFPS / FPS );
       NN = floor( length(clatvel)/Sampo);
       clatvel2 = zeros(1,NN);
       hlatvel2 = zeros(1,NN);
       for tt = 1:NN
           ta = 1+((tt-1)*Sampo);
           tb = (tt*Sampo);
           clatvel2(tt) = median( clatvel(ta:tb) );  % average, indep sets
           hlatvel2(tt) = median( hlatvel(ta:tb) );
       end
       model.latvel.cricket{tr} = clatvel2';  % down-sampled
       model.latvel.hand{tr} = hlatvel2';
    end
    TOTclatvel = [TOTclatvel ; [clatvel2' zeros(size(clatvel2')) (1:NN)']];
    TOThlatvel = [TOThlatvel ; [hlatvel2' zeros(size(hlatvel2')) (1:NN)']];
end
uuc = nanmean(TOTclatvel(:,1));
uuh = nanmean(TOThlatvel(:,1));

%% ***** look at accelerations over time
if (SGO)
   %***** located the top FP% of acceleration events
   MAXACC = 10000;
   absacc = sort( abs( TOTclatacc(:,1)) );
   Nab = length(absacc);
   thresh = absacc( floor((1-FP)*Nab) );
   %***** find high and low acc events
   Events = cell(1,2);
   VelTraces = cell(2,2);
   NL = length(TOTclatacc(:,1));
   tt = 1;
   while  (tt < NL )
       if (abs(TOTclatacc(tt,1)) > thresh)  % possible event
            tta = max(tt-FPWIN,1);
            ttb = min(tt+FPWIN,NL);
            ttz = tta:ttb;
            mpeak = max( abs(TOTclatacc(ttz,1)) );
            zpeak = find( abs(TOTclatacc(ttz,1)) == mpeak);
            if ~isempty(zpeak)
              tpeak = ttz(zpeak(1)); % time of peak
              tfr = TOTclatacc(tpeak,2);  % trial number
              tstep = TOTclatacc(tpeak,3); % frame where acc occurs
              %*** make sure time frame is within same trial
              tb = tpeak+FPWINB;
              if (tb <= NL) % && (TOTclatacc(tb,2) == tfr)  % same trial vel
                ta = tpeak-FPWINA;
                %******** categorize event 
                if ( abs(TOTclatacc(tpeak,1)) < MAXACC) && (ta >= 1)
                  %**** check if there is nan in hand within early track
                  too_early = sum(isnan(TOThlatvel(tpeak:(tpeak+floor(FPWINB/2)),1)));
                  if (1) % (too_early < (0.75*FPWINB))
                    if (TOTclatacc(tpeak,1) > 0) 
                      Events{1} = [Events{1} ; [tpeak,tfr,tstep]];
                      VelTraces{1,1} = [VelTraces{1,1} ; TOTclatvel(ta:tb,1)'];
                      VelTraces{1,2} = [VelTraces{1,2} ; TOThlatvel(ta:tb,1)'];
                    else
                      Events{2} = [Events{2} ; [tpeak,tfr,tstep]];                
                      VelTraces{2,1} = [VelTraces{2,1} ; TOTclatvel(ta:tb,1)'];
                      VelTraces{2,2} = [VelTraces{2,2} ; TOThlatvel(ta:tb,1)'];
                    end
                  end
                end
              end
              %*************
              tt = ttb;  % skip to next one
            end
       end
       tt = tt + 1;
   end
end

%%

for tr = 1:LT
    clatvel2 = model.latvel.cricket{tr}';  % down-sampled
    hlatvel2 = model.latvel.hand{tr}';
    if (1)
        %******* construct BigX matrix for regression
        LN = length(hlatvel2);
        zz = find( ~isnan(hlatvel2));
        if ~isempty(zz)
            tstart = zz(1)+1;
            % for tt = tstart:floor(BigLag/2):(LN+AntLag)
            for tt = tstart:1:(LN+AntLag)
                if (withautoreg)
                    xit = [clatvel2((tt-AntLag):-1:(tt-BigLag)) 1 hlatvel2(tt-1)];
                else
                    xit = [clatvel2((tt-AntLag):-1:(tt-BigLag)) 1 ];
                end

                yit = hlatvel2(tt);  % fit acceleration
                if ~isnan(yit)
                    BigX = [BigX ; xit];
                    BigY = [BigY ; yit];
                end
            end
        end
        %**********
    end
    %*****
    N = length(model.latvel.cricket{tr});
    model.latvel.time{tr} = (1000/FPS) * (1:N);
    %****** get the cross-correlations for included trials

    %*** zero mean, prep for correlations
    latcrick = model.latvel.cricket{tr} - uuc;
    lathand = model.latvel.hand{tr} - uuh;
    %***
    xc = nanxcov(latcrick,latcrick);
    xh = nanxcov(lathand,lathand);
    xch = nanxcov(lathand,latcrick);  % what is the right order here??
    if (length(xc) < MN)
        disp(sprintf('Throwing out short trial %d (%d frames, need %d)',...
            tr,length(xc),MN));
        LimitedMotionSet = [LimitedMotionSet, tr];  % add to throw list
    else
        mid = 1+floor(length(xc)/2);
        ctr = [(mid-MinC2):(mid+MinC2)];  % grab fix length in center
        HandCross = [HandCross ; xh(ctr)];
        CrickCross = [CrickCross ; xc(ctr)];
        HandCrick = [HandCrick ; xch(ctr)];

        %******* try to compute ARYULE filter per trace and average
        if (1)  % would need to update this with NaNs in mind
            hzz = find( ~isnan(lathand));
            czz = find( ~isnan(latcrick) );
            if (length(hzz) > (2*LAGC)) && (length(czz) > (2*LAGC))
                [arcoff,~,hrc] = aryule(lathand(hzz),LAGC);
                HandYule = [HandYule ; hrc'];
                [arcoff,~,crc] = aryule(latcrick(czz),LAGC);
                CrickYule = [CrickYule ; crc'];
            end
        end
        %**********
    end

    %******
end
%****** compute the cross correlations *****
uuhand = nanmean(HandCross);
suhand = nanstd(HandCross)/sqrt(size(HandCross,1));
uucrick = nanmean(CrickCross);
sucrick = nanstd(CrickCross)/sqrt(size(CrickCross,1));
uucross = nanmean(HandCrick);
sucross = nanstd(HandCrick)/sqrt(size(HandCrick,1));
uuhandyule = nanmean(HandYule);
suhandyule = nanstd(HandYule)/sqrt(size(HandYule,1));
uucrickyule = nanmean(CrickYule);
sucrickyule = nanstd(CrickYule)/sqrt(size(CrickYule,1));
%*********

if (1)
    %***** now go back, pre-filter traces to remove auto-correlations
    %***** and then repeat the same analyses second time
    HandCross2 = [];
    CrickCross2 = [];
    HandCrick2 = [];
    for tr = 1:LT
        %****** get the cross-correlations for included trials

        %*** zero mean, prep for correlations
        latcrick = model.latvel.cricket{tr};
        lathand = model.latvel.hand{tr};
        %*****
        patcrick = sublinearpred(latcrick,uucrickyule,0);
        pathand = sublinearpred(lathand,uuhandyule,PlotIt);
        model.patvel.cricket{tr} = patcrick;
        model.patvel.hand{tr} = pathand;
        %***
        xc = nanxcov(patcrick,patcrick);
        xh = nanxcov(pathand,pathand);
        xch = nanxcov(pathand,patcrick);  % what is the right order here??
        if (length(xc) >= MN)
            mid = 1+floor(length(xc)/2);
            ctr = [(mid-MinC2):(mid+MinC2)];  % grab fix length in center
            HandCross2 = [HandCross2 ; xh(ctr)];
            CrickCross2 = [CrickCross2 ; xc(ctr)];
            HandCrick2 = [HandCrick2 ; xch(ctr)];
        end

        %******
    end
    %***********
    %****** compute the cross correlations *****
    uuhand2 = nanmean(HandCross2);
    suhand2 = nanstd(HandCross2)/sqrt(size(HandCross2,1));
    uucrick2 = nanmean(CrickCross2);
    sucrick2 = nanstd(CrickCross2)/sqrt(size(CrickCross2,1));
    uucross2 = nanmean(HandCrick2);
    sucross2 = nanstd(HandCrick2)/sqrt(size(HandCrick2,1));
    %*********
end

%%
  
   %******* plot mean of traces triggered on ACC events
   %****** and then below that, plot the STA fits
   if (SGO > 0)
       MD = 1;  % median filter around 3 adjacent samples (+/- 1)
       hf = figure;
       % set(hf,'Position',[100 100 600 900]);
       set(hf,'Position',[100 50 600 900]);
       cvtbo = medofilt(VelTraces{1,1},MD);  % positive velocity increase only traces
       hvtbo = medofilt(VelTraces{1,2},MD);  % positive hand velocity increase only traces
       xx = (-FPWINA:FPWINB)*(1000/FPS);
       LN = size(cvtbo,1); 
       if (0)
           %***** three plots, raw velocity traces (smoothed) % positive only
           %*****              mean traces (negative reflected)
           %*****              sta from linear fit
           subplot('Position',[0.2 0.74 0.7 0.22]);
           rr = randperm(LN);
           subN = floor(1.0*LN);
           for kp = 1:subN
              rkp = rr(kp);
              plot(xx,cvtbo(rkp,:),'k-','Color',[0.7,0.7,0.7]); hold on;
              plot(xx,hvtbo(rkp,:),'k-','Color',[1.0,0.5,0.5]); hold on;
           end
           uc = nanmedian(cvtbo);
           uh = nanmedian(hvtbo);
           plot(xx,uc,'k-','LineWidth',2); hold on;
           plot(xx,uh,'r-','LineWidth',2); hold on;
           %**
           xlabel('Time (ms)');
           ylabel('Lateral Velocity (mm/s)');
           axis tight;
           V = axis;
           maxo = max(abs(V(3)),abs(V(4)));
           axis([V(1) V(2) -(maxo/2) maxo]);
           plot([V(1),V(2)],[0,0],'k-');
           plot([0,0],[-maxo,maxo],'k-');
           set(gca,'Xtick',[-50:50:150]);
           set(gca,'Ytick',[-20:20:40]);
           set(gca,'Fontsize',14);
           title('Positive Velocity Jumps');
           text(100,-0.35*maxo,sprintf('N=%d',LN),'Fontsize',14);
           text(-50,0.875*maxo,'Cricket','Fontsize',14,'Color',[0,0,0]);
           text(-50,0.700*maxo,'Hand','Fontsize',14,'Color',[1,0,0]);
       end
       
       %***************************************
       %************ then show average with negative reflected
       subplot('Position',[0.2 0.60 0.7 0.30]);
       cvtbo = medofilt([VelTraces{1,1} ; -VelTraces{2,1}],MD);  % positive velocity increase only traces
       hvtbo = medofilt([VelTraces{1,2} ; -VelTraces{2,2}],MD);  % positive hand velocity increase only traces
       uc = nanmean(cvtbo);
       nc = sum(~isnan(cvtbo));
       sc = nanstd(cvtbo) ./ sqrt(1+nc);     
       uh = nanmean(hvtbo);
       nh = sum(~isnan(hvtbo));
       sh = nanstd(hvtbo) ./ sqrt(1+nh);     
       %********
       aa = [xx fliplr(xx)];
       bb = [(uc+(nSEM*sc)) fliplr(uc-(nSEM*sc))];
       fill(aa,bb,[0,0,0],'FaceAlpha',0.3,'LineStyle','none'); hold on;
       plot(xx,uc,'k-','LineWidth',2);
       %********
       bb = [(uh+(nSEM*sh)) fliplr(uh-(nSEM*sh))];
       fill(aa,bb,[0.5,0.5,1],'FaceAlpha',0.3,'LineStyle','none'); hold on;
       plot(xx,uh,'b-','LineWidth',2,'Color',[0.3,0.3,0.6]);
       %**
       xlabel('Time (ms)');
       ylabel('Lateral Velocity (mm/s)');
       axis tight;
       V = axis;
       maxo = 25; %max(abs(V(3)),abs(V(4)));
       axis([V(1) V(2) -(maxo/2) maxo]);
       plot([V(1),V(2)],[0,0],'k-');
       plot([0,0],[-maxo,maxo],'k-');
       set(gca,'Xtick',[-50:50:150]);
       set(gca,'Ytick',[-10:10:20]);
       set(gca,'Fontsize',14);
       title('Velocity Jumps (Negative Reflected)');
       LN = size(cvtbo,1);
       text(100,-0.2*maxo,sprintf('N=%d',LN),'Fontsize',14);
       text(100,-0.35*maxo,[num2str(nSEM) ' SEM'],'Fontsize',14);
       text(-60,0.85*maxo,'Cricket','Fontsize',16,'Color',[0,0,0]);
       text(-60,0.70*maxo,'Hand','Fontsize',16,'Color',[0.3,0.3,0.6]); 
       grid off;
       %***********
       
       %********* Then show the linear prediction kernel (smooth same)
       % subplot('Position',[0.35 0.075 0.55 0.2]);
       subplot('Position',[0.20 0.10 0.70 0.30]);
       axx = (AntLag:BigLag) * (1000/FPS);
       [BigX,xx] = downsamplehist(BigX,axx,BMD);
       [b2,b2int,r2,r2int,stats] = regress(BigY,BigX);
       if (JACKN > 1)  % Jack-nife the STA (smooth it, real variance)
         bb = [];
         LN = size(BigX,1);
         for jk = 1:JACKN
           ja = 1+((jk-1)*floor(LN/JACKN));
           jb = jk*floor(LN/JACKN);
           jset = [1:ja jb:LN];
           
           b3 = regress(BigY(jset),BigX(jset,:))';
          
           bb = [bb ; b3];
         end
         b2 = mean(bb);
         b2s = std(bb) * sqrt(JACKN-1);
         b2int = [(b2+(nSEM*b2s)) ; (b2-(nSEM*b2s))];
         b2 = b2';
         b2int = b2int';
       end
       %******* get p-value on coefficients*******
       stats = regstats(BigY,BigX,'linear',{'beta','tstat'});
       BN = length(stats.beta);
       bset = [stats.tstat(1).beta(2:(BN-1)) stats.tstat(1).se(2:(BN-1)) stats.tstat(1).pval(2:(BN-1))];
       
       %***********
       aa = [xx fliplr(xx)];
       zz = 1:(size(BigX,2)-1);
       bb = [ (b2int(zz,1)') fliplr(b2int(zz,2)')];
       fill(aa,bb,[0.5,0.5,1],'FaceAlpha',0.3,'LineStyle','none'); hold on;
       plot(xx,b2(zz),'ko-','LineWidth',2); hold on;
       plot(xx,zeros(size(xx)),'k-');
       axis tight;
       V = axis;
       plot([0,0],[-0.25,0.45],'k-');
       xlabel('Time Lag (ms)');
       ylabel('Linear Coefficient');
       title('Prediction of hand from cricket velocity');
       set(gca,'Xtick',[0:50:150]);
       set(gca,'Ytick',[-0.20:0.20:0.40]);
       set(gca,'Fontsize',14);
       axis tight;
       text(100,0.35,[num2str(nSEM) ' SEM'],'Fontsize',14);
       %*******
       
   end  


return;

function [NewBigX,xx] = downsamplehist(BigX,axx,BMD)
    NewBigX = BigX;
    xx = axx;
    if (BMD <= 1)
        return;
    end
    N = size(BigX,2);
    M = floor((N-1)/BMD);
    NewBigX = zeros(size(BigX,1),(M+1));
    for tr = 1:size(BigX,1)
       nx = [];
       if (tr == 1)
          xx = [];
       end
       for j = 1:M
           na = 1+((j-1)*BMD);
           nb = j*BMD;
           nx = [nx mean(BigX(tr,na:nb))];
           if (tr == 1)
               xx = [xx mean(axx(na:nb))];
           end
       end
       NewBigX(tr,1:M) = nx;
       NewBigX(tr,(M+1)) = BigX(tr,end);
    end
return;

function ymd = medofilt(y,MD)
    ymd = y;
    N = size(y,1);
    T = size(y,2);
    wid = ceil(2*MD);
    for k = 1:N
        for t = 1:T
           ta = max(t-wid,1);
           tb = min(t+wid,T);
           % ymd(k,t) = nanmedian(y(k,ta:tb));
           ymd(k,t) = nanmean(y(k,ta:tb));
           if (0)  % Gauss smooth
             nomo = 0;
             somo = 0;
             for tt = ta:tb
               val = exp(-0.5*((tt-t)/MD)^2);
               somo = somo + (val*y(k,tt));
               nomo = nomo + val;
             end
             ymd(k,t) = (somo/nomo);
           end
        end
    end
return;

function patcrick = sublinearpred(latcrick,refcoff,plotit) 

   if (0)
     ypred = conv(refcoff,latcrick);
     LAGC = length(refcoff);
     patcrick = -ypred(1:(end-LAGC+1));
   else
     patcrick = zeros(size(latcrick));
     LAGC = length(refcoff);
     for tt = 2:length(latcrick)
         pd = 0;
         for zt = 1:LAGC
             if ((tt-zt) > 1) & ~isnan(latcrick(tt-zt))
                pd = pd + (refcoff(zt) * latcrick(tt-zt));
             end
         end
         patcrick(tt) = -pd;  % sub-tract linear prediction
      end
   end
   
   if (plotit)  % turn this on to confirm the AR fit is good
     figure(10); hold off;
     plot(latcrick,'k.-'); hold on;
     plot(patcrick,'r.-');
     input('check');
   end
   
   patcrick = latcrick - patcrick;
   
return;

%******************
function xco = nanxcov(xx,yy)
   N = length(xx);
   xco = zeros(1,((2*N)-1));
   xct = zeros(1,((2*N)-1));
   if (0)
     xx = xx - nanmean(xx);
     yy = yy - nanmean(yy);
   end
   for i = 1:N
       for j = 1:N
           tau = (i-j) + N;
           if ~isnan(xx(i)) && ~isnan(yy(j))
               xco(tau) = xco(tau) + (xx(i)*yy(j));
               xct(tau) = xct(tau) + 1;
           end
       end
   end
   xco = xco ./ xct;
return;

        