function plot_normalized_velocity_fpm(model,nq,tord)

% Generates Figure 3C from  
% Shaw,L, Wang KH, Mitchell, J (2024) Fast Prediction in Marmoset Reach-to-Grasp Movements for Dynamic Prey.
%
% inputs:
%   model = model struct of hand and cricket position
%   nq = number of interpoation points
%   tord==1:time interpolation tord==2:distance interpolation 
%
% Interpolated population average lateral and direct velocity components 
% of hands and crickets during reaches.
% Luke Shaw, Kuan Hong Wang, and Jude Mitchell 4/2024
% Matlab R2022b
%
% Reaching data structure marmo_reach_model.mat available at
% https://doi.org/10.5281/zenodo.7869286

%% define smoothing window
nTrial=size(model.x.hand,2);
sWc = 5; %smoothing window
sWh = 5;

% define some other parameters
q=1; %max hand/cricket offset val - set to 1 for no offset
for i = 1:nTrial
    Ph = [model.x.hand{i}, model.y.hand{i}];
    Pc = [model.x.cricket{i}, model.y.cricket{i}];

    Ph = smoothdata(Ph,'gaussian',sWh);
    Pc = smoothdata(Pc,'gaussian',sWc);

    nL=(size(Ph,1)-2)-q;
    k=0;

    for u=1:1:q
        k=k+1;
        nT=size(Ph,1);
        if u==1
            [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,HA_t1,CA_t1] = kinematicAnalysis(Ph,Pc);
            Cd=vecnorm(Pc(2:end,:)-Pc(1:end-1,:),2,2);
            Cdist=Cd(2:end);
            dataA{i,k} = table(VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,HA_t1,CA_t1,Cdist);
        else
            f=u-2;
            if f<nT
                Pc = [vertcat(model.x.cricketexfull{i}(end-f:end),model.x.cricket{i}(1:end-f-1)), vertcat(model.y.cricketexfull{i}(end-f:end),model.y.cricket{i}(1:end-f-1))];
                %Pc = smoothdata(Pc,'gaussian',5);
            else
                Pc = [vertcat(model.x.cricketexfull{i}(end-f:end-f+nT-1),model.x.cricket{i}(1:end-f-1)), vertcat(model.y.cricketexfull{i}(end-f:end-f+nT-1),model.y.cricket{i}(1:end-f-1))];
                %Pc = smoothdata(Pc,'gaussian',5);
            end
            Cd=vecnorm(Pc(2:end,:)-Pc(1:end-1,:),2,2);
            Cdist=Cd(2:end);
            [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,HA_t1,CA_t1] = kinematicAnalysis(Ph,Pc);
            dataA{i,k} = table(VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,HA_t1,CA_t1,Cdist);
        end

    end
end

crickdist=zeros(size(dataA,1),size(dataA,2));

for i=1:size(dataA,1)
    for j = 1:size(dataA,2)
        crickdist(i,j)=sum(dataA{i,j}.Cdist);
    end
end

%% interpolation and general velocity plot
f=1;  %future offset term

VHCa = zeros(nq,nTrial);
VHna = zeros(nq,nTrial);
VHoa = zeros(nq,nTrial);
VCoa = zeros(nq,nTrial);
HCaa = zeros(nq,nTrial);
CHa = zeros(nq,nTrial);

for i = 1:nTrial
    VHo2 = dataA{i,f}.VHo2;
    VHp2 = dataA{i,f}.VHp2;
    VHn2 = vecnorm([VHo2,VHp2],2,2); % hand velocity magnitude
    nT = length(VHo2);

    VHo1 = dataA{i,f}.VHo1;
    VHp1 = dataA{i,f}.VHp1;
    VHn1 = vecnorm([VHo1,VHp1],2,2); % hand velocity magnitude
    HA_t = dataA{i,f}.HA_t1;
    CA_t = dataA{i,f}.CA_t1;

    Dx1 = dataA{i,f}.Dx1;
    Dy1 = dataA{i,f}.Dy1;
    Dr1 = vecnorm([Dx1,Dy1],2,2); % displacement vector magnitude

    VCo1 = dataA{i,f}.VCo1; % cricket velocity othogonal to displacement vector
    VCp1 = dataA{i,f}.VCp1;

    if tord==1

        %   % by normalized time
        P = linspace(1,nT,nq);
        x = (1:nT)';
        tordname='Time';

    elseif tord==2

        %    by normalized position
        P = linspace(max(Dr1),min(Dr1),nq);
        x = Dr1;
        tordname='Dist';

    end

    VHp = interp1(x,VHp1,P);
    VHi = interp1(x,VHo1,P);
    VCi = interp1(x,VCo1,P);
    VCp = interp1(x,VCp1,P);
    VHni = interp1(x,VHn1,P); % overall hand velocity
    VHC = VHi.*sign(VCi); % hand velocity relative to cricket velocity direction
    HCi = interp1 (x,rad2deg(HA_t),P);
    CHi = interp1 (x,rad2deg(CA_t),P);
    HCa = HCi.*sign(CHi);

    VHCa(:,i) = VHC*239.76;
    VHoa(:,i) = VHi;
    VCoa(:,i) = VCi;
    VCpa(:,i) = VCp;
    VHpa(:,i) = VHp;
    VHna(:,i) = VHni;
    HCaa(:,i) = HCa;
    CHa(:,i) = CHi;
end


N=size(VHoa,2);
figure

hold on
sgtitle([tordname ' Interp: Hand and Cricket Speed']);
xp=(1:nq)';
mP = nanmean(abs(VHoa),2)*240;
sP = nanstd(abs(VHoa)*240,[],2)/sqrt(N);%./sqrt(sum(~isnan(abs(VHoa)*240),2));
del = isnan(mP);
mP(del) = [];
sP(del) = [];

p1=plot(xp,mP,'b-','linewidth',2,'DisplayName','Hand Lateral');
hold on;
fill([xp;flipud(xp)],[mP+sP;flipud(mP-sP)],'b',...
    'EdgeColor','none','FaceAlpha',0.3)
ylabel('speed (cm/s)')

mP = nanmean(abs(VHpa),2)*240;
sP = nanstd(abs(VHpa)*240,[],2)/sqrt(N);%./sqrt(sum(~isnan(abs(VHpa)*240),2));
del = isnan(mP);
mP(del) = [];
sP(del) = [];

p2=plot(xp,mP,'b--','linewidth',2,'DisplayName','Hand Direct');
hold on;
fill([xp;flipud(xp)],[mP+sP;flipud(mP-sP)],'b',...
    'EdgeColor','none','FaceAlpha',0.3)
ylabel('speed (cm/s)')

hold on

mP = nanmean(abs(VCoa),2)*240;
sP = nanstd(abs(VCoa)*240,[],2)/sqrt(N);%./sqrt(sum(~isnan(abs(VCoa)*240),2));
del = isnan(mP);
mP(del) = [];
sP(del) = [];

p3=plot(xp,mP,'r-','linewidth',2,'DisplayName','Crick Lateral');
hold on;
fill([xp;flipud(xp)],[mP+sP;flipud(mP-sP)],'r',...
    'EdgeColor','none','FaceAlpha',0.3)
ylabel('speed (cm/s)')

mP = nanmean(abs(VCpa),2)*240;
sP = nanstd(abs(VCpa)*240,[],2)/sqrt(N);%./sqrt(sum(~isnan(abs(VCpa)*240),2));
del = isnan(mP);
mP(del) = [];
sP(del) = [];

p4=plot(xp,mP,'r--','linewidth',2,'DisplayName','Crick Direct');
hold on;
fill([xp;flipud(xp)],[mP+sP;flipud(mP-sP)],'r',...
    'EdgeColor','none','FaceAlpha',0.3)
ylabel('speed (cm/s)')
ylim([0 80])
legend([p1,p2,p3,p4])


%% Functions

    function [VHp2,VHo2,VHp1,VHo1,VCp1,VCo1,Dx1,Dy1,Dxy1,HD_t1,CD_t1] = kinematicAnalysis(Ph,Pc)
        % This function extracts key kinematic parameters
        % Inputs:
        %   Ph: position of hand, [Xt,Yt], t = 1 : nT
        %   Pc: position of cricket, [Xt,Yt]

        % KH Wang 05252021

        %% number of time points
        nT = size(Ph,1);

        %% kinematic vectors
        Dch = Pc-Ph; % target displacement
        Vh = diff(Ph,1); % hand velocity
        Vc = diff(Pc,1); % cricket velocity

        % vector magnitude _radius
        Vh_r = vecnorm(Vh,2,2);
        Vc_r = vecnorm(Vc,2,2);

        % vector angles _theta
        Dch_t = atan2(Dch(:,2),Dch(:,1));
        Vh_t = atan2(Vh(:,2),Vh(:,1)); %velocity angle of hand
        Vc_t = atan2(Vc(:,2),Vc(:,1));

        % future velocity angle relative to current displacement vector
        HD_t2 = Vh_t(2:nT-1) - Dch_t(1:nT-2);

        % past velocity angle relative to current displacement vector
        HD_t1 = Vh_t(1:nT-2) - Dch_t(1:nT-2);
        CD_t1 = Vc_t(1:nT-2) - Dch_t(1:nT-2);

        %% output kinematic parameters
        % orthogonal projection of velocity vector to displacement vector
        VHo2 = sin(HD_t2).*Vh_r(2:nT-1); % future hand velocity %%
        VHo1 = sin(HD_t1).*Vh_r(1:nT-2); % past hand velocity
        VCo1 = sin(CD_t1).*Vc_r(1:nT-2); % past cricket velocity

        % parallel projection of velocity vector to displacement vector
        VHp2 = cos(HD_t2).*Vh_r(2:nT-1); % future hand velocity
        VHp1 = cos(HD_t1).*Vh_r(1:nT-2); % past hand velocity
        VCp1 = cos(CD_t1).*Vc_r(1:nT-2); % past cricket velocity

        % matched current displacement vector
        Dx1 = Dch(1:nT-2,1);
        Dy1 = Dch(1:nT-2,2);
        Dxy1= hypot(Dx1,Dy1);

    end
end
