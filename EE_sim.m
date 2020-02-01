clear all;close all
format compact

systemSize = 100;
sys = zeros(systemSize,systemSize);
occupyRate = 0.8;
occupySite0 = randsample(systemSize^2,systemSize^2*occupyRate);
year = 20000;
numt = year; 
G    = zeros(systemSize^2,numt);
C    = zeros(systemSize^2,numt);
age  = zeros(systemSize^2,numt);
biom = zeros(systemSize^2,numt);
h    = zeros(systemSize^2,numt);


J    = 1e-4;
G0   = log(0.5*J); % mean=exp(G0)=1.5e-6
sigG0= 0.9;% max=2.6e-5;min=4.1e-8;
% % r=lognrnd(G0,sigG0,[1,length(occupySite0)]);
% % histf(r,30)
hmax  = 4;
C0    = 0.4; %mean=exp(C0)=0.1
sigC0 = 0.3;%max=4.12,min=0.0003;
% % r = abs(normrnd(C0,sigC0,[1,length(occupySite0)]));
% % r=lognrnd(C0,sigC0,[1,length(occupySite0)]);
% % histf(r,30)
maxAge  = 30;
satAge  = 20;%biomass levels off with age 
seedAge = 3;%start to generate seeds
age(occupySite0,1) = datasample([10:1:maxAge],length(occupySite0));%uniform distr.
eps1  = 1.0;%effect of (C-h) on growth 
eps2  = .001; %cost of facilitation 
eps3  = 0.05; %effect of G on dh/dt
dd    = 0.005;%dispersal probably decline with distance
%greater dd, shorter dispersal
mutation = 1e-4;

maxBiom = 1;
satBiom = 0.5;%seed generation levels off with biom
k = 0.05;     %logistic formular converting age to biomass
alpha = 0.1;  %effect of h on dh/dt
beta = 0.2;   %effect of (c-h) on growth rate. 
B    = 0.05;%intrinsic growth rawte
gamma = 1.8; %what age mortality increases
seedMin = 0.2;%relationship between biomass and seed production
kk = 10;      %relationship between biomass and seed production
p_disp = 0.04^3; %threshold dispersal probability

G(occupySite0,1) = lognrnd(G0,sigG0,[1,length(occupySite0)]);
C(occupySite0,1) = lognrnd(C0,sigC0,[1,length(occupySite0)]);
biom(occupySite0,1)= maxBiom./(1+exp(-k.*(age(occupySite0,1)-satAge)));

nMutationG = 0;
nMutationC = 0;

tic
for t = 1:numt
    t
    occupySite = find(age(:,t)>0);
    G(:,t+1)    = G(:,t);
    C(:,t+1)    = C(:,t);
    
    %% seed dispersal 
    sourceLoc = occupySite(find(age(occupySite,t)>seedAge));
    sinkLoc   = find(age(:,t)==0);
    [rowSo,colSo]=ind2sub([systemSize,systemSize],sourceLoc);
    [rowSk,colSk]=ind2sub([systemSize,systemSize],sinkLoc);
    sourceXY = zeros(length(sourceLoc),2);
    sinkXY   = zeros(length(sinkLoc),2);
    sourceXY(:,1) = colSo; sourceXY(:,2) = rowSo;
    sinkXY(:,1)   = colSk; sinkXY(:,2)   = rowSk;
    distance = pdist2(sourceXY,sinkXY);
    pDisp    = 1./exp(dd.*distance);
    
    %each row is sourceLoc, each col is sinkLoc
    sourceBiom = biom(sourceLoc,t);
    pSeed = seedMin+(1-seedMin)./(1+exp(-kk.*(sourceBiom-satBiom)));
    pSeedM= repmat(pSeed,1,length(sinkLoc));
    
    sourceC    = C(sourceLoc,t);
    sinkH      = h(sinkLoc,t); 
    fitness    = repmat(sourceC,1,length(sinkLoc))-...
        repmat(sinkH',length(sourceLoc),1);
    fitness    = 1./exp(beta.*fitness.^2);
    
    prob = ((pDisp.*pSeedM).*fitness);
    [ii,jj] = max(prob); %jj is row ID of max prob
    newSites = sinkLoc(find(ii>p_disp));
    parentSite = sourceLoc(jj(find(ii>p_disp)));
    
    %mutation 
    mutationG_sites = newSites(find(mutation>rand(length(newSites),1)));
    mutationC_sites = newSites(find(mutation>rand(length(newSites),1)));
    mutatedG = lognrnd(G0,sigG0,[1,length(mutationG_sites)]);
    mutatedC = lognrnd(C0,sigC0,[1,length(mutationC_sites)]);
    nMutationG = nMutationG+length(mutationG_sites);
    nMutationC = nMutationC+length(mutationC_sites);
    
    G(newSites,t+1) = G(parentSite,t);
    G(mutationG_sites,t+1) = mutatedG;
    C(newSites,t+1) = C(parentSite,t);
    C(mutationC_sites,t+1) = mutatedC;
    biom(newSites,t)= 0.01;%initial biomass
    occupySite = unique([occupySite;newSites]);
    
    %% growth    
    rGrow = (biom(occupySite,t).*(1-biom(occupySite,t))).*...
        B.*(1 + eps1./exp(beta.*(C(occupySite,t)-h(occupySite,t)).^2)-eps2.*G(occupySite,t));
    biom(occupySite,t+1) = biom(occupySite,t) + rGrow;
    
    %% landscape change 
    h(:,t+1)=h(:,t)+(J+eps3.*(G(:,t).*biom(:,t)))./(exp(alpha*h(:,t)./max(0,hmax-h(:,t))));
    
    %% increase the age of the rest of plants 
    age(occupySite,t+1)  = age(occupySite,t)+1;
    
    %% mortality 
    m = 1./exp((maxAge-age(occupySite,t)).*gamma);
    mortLoc = find(m>rand(length(occupySite),1));
    
    biom(occupySite(mortLoc),t+1) = 0;
    G(occupySite(mortLoc),t+1)    = 0;
    C(occupySite(mortLoc),t+1)    = 0;
    age(occupySite(mortLoc),t+1)  = 0;
    
end
laps=toc


figure(3)
histogram(C(find(C(:,1)>0),1),50,'facecolor','k');hold on
histogram(C(find(C(:,1000)>0),1000),50,'facecolor','g');hold on
histogram(C(find(C(:,2000)>0),2000),50,'facecolor','b');hold on
histogram(C(find(C(:,3000)>0),3000),50,'facecolor','m');hold on
histogram(C(find(C(:,4000)>0),4000),50,'facecolor','y');hold on

% subplot(1,3,2);
histogram(G(find(G(:,1)>0),1),50,'facecolor','k');hold on
histogram(G(find(G(:,1000)>0),1000),50,'facecolor','g');hold on
histogram(G(find(G(:,2000)>0),2000),50,'facecolor','b');hold on
histogram(G(find(G(:,3000)>0),3000),50,'facecolor','m');hold on
histogram(G(find(G(:,4000)>0),4000),50,'facecolor','r');hold on

% subplot(1,3,3);
plot(sum(biom,1),'-k')
size(find(age(:,end)==0),1)/systemSize^2

figure(4)
subplot(2,2,1)
h_plot = reshape(h(:,4000),systemSize,systemSize);
surf(meshgrid(1:1:systemSize),meshgrid(1:1:systemSize)',h_plot)
subplot(2,2,2)
biom_plot = reshape(biom(:,4000),systemSize,systemSize);
surf(meshgrid(1:1:systemSize),meshgrid(1:1:systemSize)',biom_plot)
subplot(2,2,3)
G_plot = reshape(G(:,4000),systemSize,systemSize);
surf(meshgrid(1:1:systemSize),meshgrid(1:1:systemSize)',G_plot)
subplot(2,2,4)
C_plot = reshape(C(:,4000),systemSize,systemSize);
surf(meshgrid(1:1:systemSize),meshgrid(1:1:systemSize)',C_plot)


%%
size(unique(C(:,end)))
size(unique(G(:,end)))

ecoDebt= zeros(2001,1);
for k = 1:2001
    ecoDebt(k) = mean(C(:,k)-h(:,k)); 
end
plot(ecoDebt)

figure(1)
plot(mean(C,1),'-k');hold on
plot(mean(h,1),'-b');grid on

plot(G(:,end),h(:,end),'ok')
plot(mean(age,1))

plot(C(:,4000),G(:,4000),'ok')

plot(age(:,end),biom(:,end),'ok')
