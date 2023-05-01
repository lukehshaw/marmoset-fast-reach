function grasp_cluster_mfr()

% Generates Figure 2D,E,F. Clusters z-scored finger separation time series 
% and plots clustered finger separation, velocity, and principle components. 
% Requires marmo_grasp_model.mat. Marmo_grasp_model.mat contains reaches
% with acceptable finger tracking and extended tracking to capture the
% beginning of hand closure.
%
% Shaw,L, Wang KH, Mitchell, J (2023) Fast Prediction in Marmoset Reach-to-Grasp Movements for Dynamic Prey.
%
% Jude Mitchell, Kuan Hong Wang, and Luke Shaw 4/2023
% MATLAB R2022b
%
% Reaching data structure marmo_reach_model.mat available at
% https://doi.org/10.5281/zenodo.7869286


%% Load Data
load marmo_grasp_model.mat

%% Set these parameters
MC=3; %max number of cluser
tord=1; % interp method: 1 = TIME 2 = DISTANCE from first cricket loc 3=Distance from last hand loc 4=Distance from last cricket loc
VorHS=[2 3]; %plots 1=z scores 2 = velocity 3= nonz-scored hand spread traces as function of clustering
mdwtorkmean=1; % 1=mdwtcluster 2=kmean

colors=[[0.4660 0.6740 0.1880 .2];[0.4940 0.1840 0.5560 .2];[0.9290 0.6940 0.1250 .2]];


%% Separate into knuckles and tips

for r=1:size(work.x.food,2)

    for j=1:length(work.x.food{r})
        f=1;
        for e=1:4
            AB=norm([work.x.fingerk{r}(j,e),work.y.fingerk{r}(j,e)]-[work.x.fingerk{r}(j,e+1),work.y.fingerk{r}(j,e+1)]);
            work.fdistk{r}(j,e)=AB;
        end
    end
    work.fdistfullk{r}=nanmean(smoothdata(work.fdistk{r}(:,2:4),1,'gaussian',10),2);
end


for r=1:size(work.x.food,2)
    for j=1:length(work.x.food{r})
        f=1;
        for e=1:4

            AB=norm([work.x.fingert{r}(j,e),work.y.fingert{r}(j,e)]-[work.x.fingert{r}(j,e+1),work.y.fingert{r}(j,e+1)]);
            work.fdistt{r}(j,e)=AB;
        end
    end
    work.fdistfullt{r}=nanmean(smoothdata(work.fdistt{r}(:,2:4),1,'gaussian',10),2);
end

%% Interpolate and generate Z scored versions
nq=51;
int_limit=1; %set higher than 1 if you want to downsample
Dist=[];
N=size(work.x.hand,2);

for i = 1:size(work.x.hand,2)
    ints=1:int_limit:(length(work.x.food{i})-1);

    if tord==1

        %   % by normalized time
        P = linspace(1,length(work.x.food{i})-1,nq);
        x = (1:length(work.x.food{i})-1)';
        tordname='time';

    elseif tord==2

        %    by normalized position relative to cricket
        xD=work.x.hand{i}(ints)-ones(length(ints),1)*work.x.food{i}(1);
        yD=work.y.hand{i}(ints)-ones(length(ints),1)*work.y.food{i}(1);
        dis=hypot(xD,yD);
        dist{i}=fillmissing(dis,'linear');
        P = linspace(max(dist{i}),min(dist{i}),nq);
        x = dist{i};
        tordname='distcrick';

    elseif tord==3

        %    by normalized position relative to hand
        xD=work.x.hand{i}(ints)-ones(length(ints),1)*work.x.hand{i}(end);
        yD=work.y.hand{i}(ints)-ones(length(ints),1)*work.y.hand{i}(end);
        dis=hypot(xD,yD);
        dist{i}=fillmissing(dis,'linear');
        P = linspace(max(dist{i}),min(dist{i}),nq);
        x = dist{i};
        tordname='distendhand';

    elseif tord==4
        if sum(isnan(work.x.food{i}))==0
            nanCrick(i)=length(work.x.food{i})-1;
            continue
        end
        nanCrick(i)=min(find(isnan(work.x.food{i})==1))-1;

        xD=work.x.hand{i}(ints)-ones(length(ints),1)*work.x.food{i}(nanCrick(i));
        yD=work.y.hand{i}(ints)-ones(length(ints),1)*work.y.food{i}(nanCrick(i));
        dis=hypot(xD,yD);
        dist{i}=fillmissing(dis,'linear');
        P = linspace(max(dist{i}),min(dist{i}),nq);
        x = dist{i};
        tordname='distcrick';

    end

    Kix{i}=interp1(x',work.x.fingerk{i}(1:end-1,:),P,'linear','extrap');
    Kiy{i}=interp1(x',work.y.fingerk{i}(1:end-1,:),P,'linear','extrap');
   % if find(all(isnan(work.x.fingert{i}),1)==1)
        nanfind=all(isnan(work.x.fingert{i}),1);
        Tix{i}=NaN(nq,5);
        Tiy{i}=NaN(nq,5);
        Tix{i}(:,nanfind==0)=interp1(x',work.x.fingert{i}(1:end-1,nanfind==0),P,'pchip','extrap');
        Tiy{i}(:,nanfind==0)=interp1(x',work.y.fingert{i}(1:end-1,nanfind==0),P,'pchip','extrap');
%     else
%         Tix{i}=interp1(x',work.x.fingert{i}(1:end-1,:),P,'pchip','extrap');
%         Tiy{i}=interp1(x',work.y.fingert{i}(1:end-1,:),P,'pchip','extrap');
%     end
    Hix{i}=interp1(x',work.x.hand{i}(1:end-1,:),P,'linear','extrap');
    Hiy{i}=interp1(x',work.y.hand{i}(1:end-1,:),P,'linear','extrap');
    Fix{i}=interp1(x',work.x.food{i}(1:end-1,:),P,'linear');
    Fiy{i}=interp1(x',work.y.food{i}(1:end-1,:),P,'linear');

    FK(i,:)=interp1(x',fillmissing(work.fdistfullk{i}(2:end),'linear')',P);
    FKz(i,:)=zscore(interp1(x',fillmissing(work.fdistfullk{i}(2:end),'linear')',P));
    FT(i,:)=interp1(x',fillmissing(work.fdistfullt{i}(2:end),'linear')',P);
    FTz(i,:)=zscore(interp1(x',fillmissing(work.fdistfullt{i}(2:end),'linear')',P));
    K3(i,:)=interp1(x',fillmissing(work.fdistk{i}(2:end,3),'linear')',P);
    K3z(i,:)=zscore(interp1(x',fillmissing(work.fdistk{i}(2:end,3),'linear')',P));
    T3(i,:)=interp1(x',fillmissing(work.fdistt{i}(2:end,3),'linear')',P);
    T3z(i,:)=zscore(interp1(x',fillmissing(work.fdistt{i}(2:end,3),'linear')',P));
    K2(i,:)=interp1(x',fillmissing(work.fdistk{i}(2:end,2),'linear')',P);
    K2z(i,:)=zscore(interp1(x',fillmissing(work.fdistk{i}(2:end,2),'linear')',P));
    T2(i,:)=interp1(x',fillmissing(work.fdistt{i}(2:end,2),'linear')',P);
    T2z(i,:)=zscore(interp1(x',fillmissing(work.fdistt{i}(2:end,2),'linear')',P));
    VI(i,:)=interp1(x',smoothdata(work.xy.handv{i},'gaussian',5)',P);
    VIz(i,:)=zscore(interp1(x',smoothdata(work.xy.handv{i},'gaussian',5)',P));
end

mFK=nanmedian(FK);
mFKz=nanmedian(FKz);
mFT=nanmedian(FT);
mFTz=nanmedian(FTz);
mVI=nanmean(VI);
mT3=nanmean(T3);
mK3=nanmean(K3);
sFK=nanstd(FK)/sqrt(N);
sFT=nanstd(FT)/sqrt(N);
sVI=nanstd(VI)/sqrt(N);


%% Cluster Z score knuckles or tips

if mdwtorkmean==1
    cluster=mdwtcluster(FKz,'maxclust',MC);
elseif mdwtorkmean==2
    cluster.IdxCLU=kmeans(FKz,MC,'Replicates',10);
end

for s=VorHS

    figure
    t=tiledlayout(1,MC);
    for i=1:MC
        nexttile
        if s==1
            plot(FKz(cluster.IdxCLU(:,1)==i,:)','Color',colors(i,:)); hold on;
            plot(nanmean(FKz(cluster.IdxCLU(:,1)==i,:)), 'k','LineWidth',4,'Color',colors(i,1:3))
            ylim([-2 3])
            xlim([1 nq])
            set(gca,'TickLength',[0 .01])
            subtitle('z score separation')
        elseif s==2
            plot(VI(cluster.IdxCLU(:,1)==i,:)','Color',colors(i,:)); hold on;
            plot(nanmean(VI(cluster.IdxCLU(:,1)==i,:)), 'k','LineWidth',4,'Color',colors(i,1:3))
            ylim([0 90])
            xlim([1 nq])
            set(gca,'TickLength',[0 .01])
            subtitle('velocity')
            t.Title.String='Fig 2F : Cluster Group Hand Velocity';
            t.Title.FontWeight='bold';
        elseif s==3
            plot(FK(cluster.IdxCLU(:,1)==i,:)','Color',colors(i,:)); hold on;
            plot(nanmean(FK(cluster.IdxCLU(:,1)==i,:)), 'k','LineWidth',4,'Color',colors(i,1:3))
            xlim([1 nq])
            ylim([.25 .65])
            set(gca,'TickLength',[0 .01])
            subtitle('separation')
            t.Title.String='Fig 2E : Cluster Group Finger Separation';
            t.Title.FontWeight='bold';
        end
        chosen{i}=work.check(cluster.IdxCLU(:,1)==i,:);
        title(num2str(i));
        set(gcf,'color','white')
    end
end


%% PLOT ONTO PCA SPACE

[coeff,score]=pca(FKz);
figure

for i=1:MC
    zz=find(cluster.IdxCLU(:,1)==i); hold on
    scatter(score(zz,1),score(zz,2),'filled','MarkerFaceColor',colors(i,(1:3)));
    Info{i}=['cluster ' num2str(i)];
end

for n = 1:size(score,1)
    if n<10
        Num=['0' num2str(n)];
    else
        Num=num2str(n);
    end
    text(score(n,1)-.15,score(n,2)+.4,Num)
end

set(gcf,'color','white')
legend(Info,'Location','southwest')
xlabel('PC1');
ylabel('PC2');
title('Fig 2D: Finger Separation Clustering')

end