function regression_model_plot_mfr(res)

% Generates Figure 4E plot from multivariate regression analysis. Requires
% that regression_model_mfr has been run. Enter resolution value.
%
% inputs:  
%   res     Resolution of delay and prediction sampling used in
%           regression_model_mfr
%               1 = high resolution (~4ms)
%               2 = medium resolution (~12ms) 
%               3 = low resolution (~20ms)
%
% Shaw,L, Wang KH, Mitchell, J (2023) Fast Prediction in Marmoset Reach-to-Grasp Movements for Dynamic Prey.
%
% Jude Mitchell, Kuan Hong Wang, and Luke Shaw 4/2023
% MATLAB R2022b
%
% Reaching data structure marmo_reach_model.mat available at
% https://doi.org/10.5281/zenodo.7869286

%%

load(['marm_regress_' num2str(res) '.mat']); %saved output variable space from regression_model_fpm

%% 

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
set(gcf,'color','white');
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
colorbar('Ticks',[200,250,300,350,400],'TickLabels',{'-200','-250','-300','-350','-400'}); %switch back to negative values
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

end
