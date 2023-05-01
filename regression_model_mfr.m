function regression_model_mfr(res)

% Generates data to plot Figure 4E using the equation in Figure 4C.
% Saves output variables required for plotting as 'output' workspace
% Plot script is regression_model_plot_mfr.m
% Shaw,L, Wang KH, Mitchell, J (2023) Fast Prediction in Marmoset Reach-to-Grasp Movements for Dynamic Prey.

% inputs:  
%   res     Resolution of delay and prediction sampling
%               1 = high resolution (~4ms)
%               2 = medium resolution (~12ms) 
%               3 = low resolution (~20ms)
% 
% Will save a .mat file 'marm_reg_(n=res value).mat' for later plotting without running this script  
%
% Eg. regression_model_mfr(3); for low resolution
%
% This function performs multivariate linear regression modeling according
% to the equation in Figure 4C and outputs AIC difference.
% Jude Mitchell, Kuan Hong Wang, and Luke Shaw 4/2023
% MATLAB R2022b
%
% Reaching data structure marmo_reach_model.mat available at
% https://doi.org/10.5281/zenodo.7869286

%% import model data
load('marmo_reach_model.mat','model');

%% 
output = ['marm_regress_' num2str(res) '.mat'];

%%
if res==1
    q=1:37;
    p=-60:60;
else
    q=1:(3*(res-1)):37; %visuomotor delay values to test 
    p=-60:(5*(res-1)):60; %prediction tau values to test
end

%%%%%

nTrial = 1:size(model.x.cricket,2);

AICj=ones(5,length(q),length(p));
BClat=ones(length(q),length(p));
sumTbl=cell(length(q),length(p));

tic
if (1)
  disp('Running main model fit');
  [AICj,BClat,sumTbl,nT]=zkinematicRegressionAIC(model,q,p);
  min(min(squeeze(AICj(2,:,:))))
end
toc
disp('Finished the main model fit, continuing on bootstraps');

%******** Do Jacknife on dataset to get error bars on estimators
if (1)
  JackN = 4; %10; % size(model.id,1);   % trial numbers
  ThrowN = ceil( size(model.id,1)/JackN);
  AICjack = cell(1,JackN);
  for bk = 1:JackN
   disp(sprintf('Jacknife %d of %d',bk,JackN));
   %****** adjust model to include a subset
   jmodel = model;
   %****** remove a subset of trials here
   astart = 1 + ((bk-1)*ThrowN);
   afini = min(size(model.id,1),(bk*ThrowN));
   [astart,afini]
   for jk = astart:afini
       jmodel.x.hand{jk} = [];  % skip this trial of data
   end
   %******************
   tic
   [AICja,~,~,~]=zkinematicRegressionAIC(jmodel,q,p);
   toc
   AICjack{bk} = AICja;
   min(min(squeeze(AICja(2,:,:))))
  end
  disp('Finished jacknife fits');  
  save(output);
end

%% Plot

AICbase = squeeze(AICj(1,:,:) - AICj(4,:,:));  % base model, normed auto model (4)
basebest = min(min(AICbase));
AICpred = squeeze(AICj(2,:,:) - AICj(4,:,:));
AICdiff = AICpred - AICbase;
%***** works out the same as AICdiff
if (0)
  zmin = find( AICdiff == min(min(AICdiff)) );  
else
  %***** determin the min on the visuomotor delay of 80 ms
  X = q ./ 240;
  Xdist = (X - 0.08) .^ 2;
  mdist = min(Xdist);
  zz = find( Xdist == mdist);
  xmin = zz(1);
  mino = min( squeeze(AICpred(xmin,:)) );
  zmin = find( AICpred == mino );  
end
%****************
X=repmat(q./240,length(p),1)';
Y=repmat(p./240,length(q),1);
%*********

%***** store size of nT in next round
nT = 0;
for jk = 1:length(model.x.hand)
    nT = nT + length(model.x.hand{jk});
end

%****** compute Jacknifed version
for jk = 1:JackN
   AICjakbase(jk,:,:) = squeeze( (AICjack{jk}(1,:,:) - AICjack{jk}(4,:,:)) );
   AICjakpred(jk,:,:) = squeeze( (AICjack{jk}(2,:,:) - AICjack{jk}(4,:,:)) );
   % AICjakdiff(jk,:,:) = squeeze( (AICjack{jk}(2,:,:) - AICjack{jk}(1,:,:)) );  % not identical at zero
   AICjakdiff(jk,:,:) = squeeze( (AICjack{jk}(2,:,:) - AICjack{jk}(5,:,:)) ); % converges at zero 
end
AIC_base_std = squeeze( std(AICjakbase) ) * sqrt(JackN-1);
AIC_pred_std = squeeze( std(AICjakpred) ) * sqrt(JackN-1);
AIC_diff_std = squeeze( std(AICjakdiff) ) * sqrt(JackN-1);
%AICstd = AIC_diff_std;
AICstd = 0.5*sqrt(AIC_diff_std.^2+AIC_pred_std.^2); 
%This choice of error bound estimation lies between the best motivated and
%lowest std (model difference) and an over-estimate of std from the
%prediction model alone

%********* colormap to emphasize trough
cmap = viridis(101);
%***************
iX = q/240;
iY = p/240;
PeakVM = X(zmin);
PeakTP = Y(zmin);  
xp = find( PeakVM == iX );
yp = find( PeakTP == iY );
%****** modify to show prediction performance
TPcurve = -AICpred(xp,:);
TPstd = AICstd(xp,:);  % shows variance pred model to PNP line
TPpnp = mean(-AICbase(xp,:));

%***********
VMcurve = AICdiff(:,yp)';
VMstd = AICstd(:,yp)';  % from 10 fold Jacknife
%******
%***** model with no cricket term, what is best prediction?
VMDcurve = AICbase(:,yp)';
VMDstd = AIC_base_std(:,yp)';
%*********

hf = figure;
set(hf,'Position',[100 100 800 800]);
subplot('Position',[0.35 0.10 0.5 0.8]);
diffX = 0.5*median(diff(iX));
diffY = 0.5*median(diff(iY));
minx = min(iX)-diffX;
maxx = max(iX)+diffX;
miny = min(iY)-diffY;
maxy = max(iY)+diffY;
imagesc(flipud(iX),iY,flipud(-AICpred')); hold on;
plot([minx,maxx],[0,0],'k:','LineWidth',1.5);
plot([PeakVM,PeakVM],[miny,maxy],'k--','LineWidth',2);
plot([minx,maxx],[-PeakTP,-PeakTP],'k--','LineWidth',2);
%** box in image
plot([minx,minx,maxx,maxx,minx],[miny,maxy,maxy,miny,miny],'k-',...
        'LineWidth',1.5);
text((minx+(1.275*(maxx-minx))),(miny+(0.3*(maxy-miny))),...
        'Goodness of Fit','Rotation',-90,'Fontsize',16);
text((minx+(1.275*(maxx-minx))),(miny+(0.0*(maxy-miny))),...
        'Better','Fontsize',16);
text((minx+(1.275*(maxx-minx))),(miny+(1.0*(maxy-miny))),...
        'Worse','Fontsize',16);
axis off;
%*** bound image with a line
title('Figure 4E: Predictive Model');
set(gca,'Fontsize',14);
%******** make a special color map here
colormap(cmap);
colorbar;
%****** add horizontal axis of surface plot
subplot('Position',[0.35 0.09 0.425 0.001]);
axis([min(iX) max(iX) 0 1]);
xticks([-0.04 0 0.04 0.08 0.12]);
xticklabels({'-40','0','40','80','120'});
yticks([]);
set(gca,'Fontsize',14);
set(gca,'Linewidth',2);
xlabel('Visuomotor Delay (secs)')
%********    

%******* Plot marginal on the prediction term
subplot('Position',[0.15 0.10 0.190 0.8]);
SEM = 1;
aminy = min(TPcurve);
amaxy = max(TPcurve);
adiff = (amaxy-aminy);
aminy = aminy - 0.30*adiff;
amaxy = amaxy + 0.35*adiff;
fill([aminy,aminy,amaxy,amaxy,aminy],[miny,maxy,maxy,miny,miny],[1,1,1]); hold on;
plot([aminy,aminy,amaxy,amaxy,aminy],[min(iY),max(iY),max(iY),min(iY),min(iY)],'k-','LineWidth',2); hold on;
%******
xaa = [iY fliplr(iY)];
yaa = [(TPcurve+(SEM*TPstd)) fliplr(TPcurve-(SEM*TPstd))];
fill(yaa,xaa,[0.5,0.5,0.5],'Linestyle','none','FaceAlpha',0.5);
plot(TPcurve,iY,'k-','LineWidth',2); hold on;
%************
PP2Dcolo = [0,0,0];
PN2colo = [0.8,0.0,0.4];
plot([340,340],[miny,maxy],'k:','LineWidth',2,'Color',PN2colo);

plot([aminy,amaxy],[0,0],'k:','LineWidth',1.5);
plot([aminy,amaxy],[PeakTP,PeakTP],'k--','LineWidth',2);
axis([aminy amaxy min(iY) max(iY)]);
%*** make your own xticks above
xticks([]);
text(aminy+(0.21*(amaxy-aminy)),miny+(1.08*(maxy-miny)),'Goodness','Fontsize',16);
text(aminy+(0.34*(amaxy-aminy)),miny+(1.03*(maxy-miny)),'of Fit','Fontsize',16);
%*******
text(aminy-(0.60*(amaxy-aminy)),miny+(1.00*(maxy-miny)),'Future','Fontsize',16);
text(aminy-(0.60*(amaxy-aminy)),miny+(0.00*(maxy-miny)),'Past','Fontsize',16);
%*************
yticks([-0.2 -0.1 0.0 0.1 0.2]);
yticklabels({'-200','-100','0','100','200'}); 
set(gca,'Fontsize',14);
ylabel('Prediction Tau (ms)');
set(gcf,'color','white');

%%
    function [AICj,BClat,sumTbl,nT] = zkinematicRegressionAIC(model,q,p)

        nTrial = 1:size(model.x.cricket,2);


        for k=1:length(q) %delay term
            for j=1:length(p) %predict term
                %*********
                %****** doing this an older way without table to make sure it works
                zVHp2 = [];
                zVHo2 = [];
                zVHp1 = [];
                zVHo1 = [];
                zVCp1 = [];
                zVCo1 = [];
                zDx1 = [];
                zDy1 = [];
                zDxy1 = [];
                zPNx = [];
                zPNy = [];
                zPPx = [];
                zPPy = [];
                %********

                for g = 1:length(nTrial)

                    if isempty(model.x.hand{g})
                        continue;  % used in Jacknifing to skip subsets of trials
                    end

                    Ph = [model.x.hand{g}, model.y.hand{g}];
                    Pc = [model.x.cricket{g}, model.y.cricket{g}];
                    Ph = smoothdata(Ph,'gaussian',5);
                    Pc = smoothdata(Pc,'gaussian',5);
                    nT=size(Ph,1);

                    if q(k)==1
                        [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,PNx,PNy,PPx,PPy] = zkinematicAnalysis(Ph,Pc);
                    else
                        if (q(k) >= 2)
                            f=q(k)-2; %pulls from past cricket data (exfull) and aligns with hand data
                            if f < nT %if vm delay less than length of reach
                                Pc = [vertcat(model.x.cricketexfull{g}(end-f:end),model.x.cricket{g}(1:end-f-1)), ...
                                    vertcat(model.y.cricketexfull{g}(end-f:end),model.y.cricket{g}(1:end-f-1))];
                                Pc = smoothdata(Pc,'gaussian',5);
                                % may be an issue here, Ph does not shift back for
                                % computing the past range vector, Dx1 and Dy1
                            else
                                Pc = [vertcat(model.x.cricketexfull{g}(end-f:end-f+nT-1),model.x.cricket{g}(1:end-f-1)), ...
                                    vertcat(model.y.cricketexfull{g}(end-f:end-f+nT-1),model.y.cricket{g}(1:end-f-1))];
                                Pc = smoothdata(Pc,'gaussian',5);
                            end
                            [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,PNx,PNy,PPx,PPy] = zkinematicAnalysis(Ph,Pc);
                        else
                            f = -q(k)+1;
                            if f < (nT-1) %if vm delay less than length of reach
                                Ph = [model.x.hand{g}(1:(end-f)), model.y.hand{g}(1:(end-f))]; % hand steps into past
                                Ph = smoothdata(Ph,'gaussian',5);
                                Pc = [model.x.cricket{g}((f+1):end), model.y.cricket{g}((f+1):end)];  % cricket moves into future
                                Pc = smoothdata(Pc,'gaussian',5);
                            else
                                continue;  % skip this trial then, not enough data
                            end
                            [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,PNx,PNy,PPx,PPy] = zkinematicAnalysis(Ph,Pc);
                        end
                    end

                    VCo1=VCo1*p(j);   % is VCo1 in (mm/s) ... what about step size?
                    VCp1=VCp1*p(j);


                    %*** only include beyond some TM from start of trial
                    TM = 0;
                    if (TM)
                        leno = length(VHp2);
                        if leno > TM
                            xx = TM:leno;
                            VHp2 = VHp2(xx,1);
                            VHo2 = VHo2(xx,1);
                            VHp1 = VHp1(xx,1);
                            VHo1 = VHo1(xx,1);
                            VCp1 = VCp1(xx,1);
                            VCo1 = VCo1(xx,1);
                            Dx1 = Dx1(xx,1);
                            Dy1 = Dy1(xx,1);
                            Dxy1 = Dxy1(xx,1);
                            PNx = PNx(xx,1);
                            PNy = PNy(xx,1);
                            PPx = PPx(xx,1);
                            PPy = PPy(xx,1);
                        else
                            VHp2 = [];
                            VHo2 = [];
                            VHp1 = [];
                            VHo1 = [];
                            VCp1 = [];
                            VCo1 = [];
                            Dx1 = [];
                            Dy1 = [];
                            Dxy1 = [];
                            PNx = [];
                            PNy = [];
                            PPx = [];
                            PPy = [];
                        end
                    end

                    %********
                    zVHp2 = [zVHp2; VHp2];
                    zVHo2 = [zVHo2; VHo2];
                    zVHp1 = [zVHp1; VHp1];
                    zVHo1 = [zVHo1; VHo1];
                    zVCp1 = [zVCp1; VCp1];
                    zVCo1 = [zVCo1; VCo1];
                    zDx1 = [zDx1; Dx1];
                    zDy1 = [zDy1; Dy1];
                    zDxy1 = [zDxy1; Dxy1];
                    zPNx = [zPNx ; PNx];
                    zPNy = [zPNy ; PNy];
                    zPPx = [zPPx ; PPx];
                    zPPy = [zPPy ; PPy];
                    %******

                    %******* Don't understand how this is working ...
                    %  but appears that T{k+1} concats all T{k} before it
                    T{g} = table(VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,PNx,PNy,PPx,PPy);
                    %*****
                end

                %**********
                VHp2 = zVHp2;
                VHo2 = zVHo2;
                VHp1 = zVHp1;
                VHo1 = zVHo1;
                VCp1 = zVCp1;
                VCo1 = zVCo1;
                Dx1 = zDx1;
                Dy1 = zDy1;
                Dxy1 = zDxy1;
                PNx = zPNx;
                PNy = zPNy;
                PPx = zPPx;
                PPy = zPPy;
                BB = table(VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,PNx,PNy,PPx,PPy);
                %************

                % sanity check, this code gives the same results, but will be
                % robust to the blanks on the Jacknife subsets
                dTB = cat(1,BB);
                [AICj(:,k,j),BClat(k,j),sumTbl{j,k},nT(j,k)]=zkinematicRegressionInHandDxyRJack(dTB);
            end
            disp(sprintf('Testing delay term %d of %d with nT = %d',q(k),q(end),nT(end,k)));
        end
    end

%%

function [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,PNx,PNy,PPx,PPy] = zkinematicAnalysis(Ph,Pc)

%   Ph: position of hand, [Xt,Yt], t = 1 : nT
%   Pc: position of cricket, [Xt,Yt]

% KH Wang 05252021 + LH Shaw edits

% number of time points

nT = size(Ph,1);

% kinematic vectors
Dch = Pc-Ph; % target displacement (note Ph is not time shifted here!)
Vh = diff(Ph,1); % hand velocity
Vc = diff(Pc,1); % cricket velocity

%***** does velocity need to be in (mm/sec)?  if so, times step size?
if (0)   % this would be essential to the steering models (right coordinates)
         % but possibly not to the 2D models where mvregress can compensate
  Vh = Vh * 240;
  Vc = Vc * 240;
end

%*********************************** JUDE ADDITION FOR PNP MODEL
% proportinal navigation commands
for k = 2:(length(Dch)-1)
   r = [Dch(k,:) 0];  % range vector, go 3D to keep math same with Brighton&Taylor
   r2 = [Dch(k-1,:) 0];
   dr = r-r2;
   vc = [Vc(k-1,:) 0];
   vh = [Vh(k-1,:) 0];
   %****
   lambdaDot = cross(r2 , (vc-vh))/sum(r2.^2);
   %****** other 2D form of prop nav
   %Calculates commanded acceleration under delayed PN guidance law
   birdAcc2 = cross(lambdaDot,vh);
   
   % it include this, equivalent to mixed strategy, becomes prediction
   % model when fit with mvregress because parametes are unconstrained
   % birdAcc2 = dr; % lambdaDot * dr;
   PNx(k-1,1) = birdAcc2(1);
   PNy(k-1,1) = birdAcc2(2);
   
   %*********
   delta=acos(dot(r2,vh)/(norm(r2)*norm(vh)));
    
   %Creates unit vector normal to line of sight (r) and velocity (v)
   n=cross(r2,vh)/(norm(r2)*norm(vh)); 
   n=n/norm(n); 
    
   %Calculates commanded acceleration under delayed PP guidance law
   birdAcc = delta*cross(n,vh);  % command for PP
   
   PPx(k-1,1) = birdAcc(1);
   PPy(k-1,1) = birdAcc(2);
   
end
%**************************************

% output kinematic parameters
% orthogonal projection of velocity vector to displacement vector
VHo2 = Vh(2:nT-1,1); % future hand velocity %%
VHo1 = Vh(1:nT-2,1); % past hand velocity
VCo1 = Vc(1:nT-2,1); % past cricket velocity

% parallel projection of velocity vector to displacement vector
VHp2 = Vh(2:nT-1,2); % future hand velocity
VHp1 = Vh(1:nT-2,2); % past hand velocity
VCp1 = Vc(1:nT-2,2); % past cricket velocity

% matched current displacement vector 
Dx1 = Dch(1:nT-2,1);
Dy1 = Dch(1:nT-2,2);
Dxy1= hypot(Dx1,Dy1);

end


%%
function [AICstats,BoC,sumTbl,nT] = zkinematicRegressionInHandDxyRJack(dT)
% regression modeling on kinematic features
% displacement scalar (dxy1) rather than dx1 and dy1

%   dT: table containing predictors and response variables for regression

% regression and hypothesis testing 
nT = size(dT,1);
     
% response variable
Y = [dT.VHp2, dT.VHo2];
% predictor varialbe. Two models
X = cell(2,1); 
X{1} = [dT.VHp1, dT.VHo1]; % without cricket
X{2} = [[dT.Dy1, dT.Dx1] + [dT.VCp1, dT.VCo1], X{1}]; % with cricket

% variable names
varNames = [dT.Properties.VariableNames,'Cnst'];
varNamesP = cellfun(@(x) {[x,'_p']},varNames);
varNamesO = cellfun(@(x) {[x,'_o']},varNames);
predNames = cell(4,1);
predNames{1} = [varNamesP([3 4]), varNamesO([3 4])]; %without cricket
predNames{2} = [varNamesP([5 6 3 4]), varNamesO([5 6 3 4])]; %with cricket

%** replace X{1} so it contains range vector as well
if (1)
  X{3} = X{1};  % momentum term only, use to normalize other errors
  predNames{3} = [varNamesP([3 4]), varNamesO([3 4])]; %without cricket
  %****
  %*** best 2D model with no prediction term
  X{4} = [[dT.Dy1, dT.Dx1], X{1}];  % same, but no cricket => tau_pred = 0
  predNames{4} = [varNamesP([5 6 3 4]), varNamesO([5 6 3 4])]; %with cricket
  %*** only PN model (reactive)
  % X{1} = [[dT.PNy, dT.PNx], X{1}]; 
  % predNames{1} = [varNamesP([5 6 3 4]), varNamesO([5 6 3 4])]; %with cricket
  %***********
  %***  mixed PP and PN model (best reactive)
  % X{1} = [[dT.Dy1, dT.Dx1],[dT.PNy, dT.PNx], X{1}]; 
  X{1} = [[dT.PPy, dT.PPx],[dT.PNy, dT.PNx], X{1}]; 
  predNames{1} = [varNamesP([5 6 3 4 1 2]), varNamesO([5 6 3 4 1 2])]; %with cricket
 
end


for m = 1:4  % m = 1 does not include cricket x,y (range only), m = 2 has both
    
    Xd = X{m}; 
    
    [B,~,~,CovB,logL] = mvregress(Xd,Y); %%% MVREGRESS IS HERE!!
    
    %***** the neg log likelihood will depend on N, normalize that
    % logL = logL/length(Y);  % normalize by number of observations included
    %***************
    
    PartialStdErr = diag(sqrt(CovB));
    tRatio = B(:)./PartialStdErr;
    pVals = 2*(1-tcdf(abs(tRatio),nT-2));
    pVals = reshape(pVals,[],2);% two columns for XY 
    AIC(m) = 2*numel(B)-2*logL;
    %AIC(m) = logL;

    varN = reshape(predNames{m},[],2);
    sumTbl = table(varN,B,pVals,'VariableNames',{'predictor','beta','pVal'});
    BoC=sumTbl.beta(2,2);
end

AICstats=[AIC(1) AIC(2) (AIC(2)-AIC(3)) AIC(3) AIC(4)];

end

end
