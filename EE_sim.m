clear all;close all

systemSize = 100;
sys = zeros(systemSize,systemSize);
occupyRate = 0.8;
occupySite0 = randsample(systemSize^2,systemSize^2*occupyRate);

year = 2000;
dt = 1;
numt = year; 
G0   = -7; %mean=exp(G0)=0.0010
sigG0= 0.5;%max=0.0076;min=0.00012;
r=lognrnd(G0,sigG0,[1,length(occupySite0)]);
hmax = 4;
C0   = -3; %mean=exp(C0)=0.05
sigC0= 1.1;%max=0.51,min=0.005;

maxAge  = 200;
satAge  = 70;%biomass levels off with age 
seedAge = 10;%start to generate seeds
maxBiom = 1;
satBiom = 0.5;%seed generation levels off with biom
k = 0.05;     %logistic formular converting age to biomass


G    = zeros(systemSize^2,numt);
C    = zeros(systemSize^2,numt);
age  = zeros(systemSize^2,numt);
biom = zeros(systemSize^2,numt);
h    = zeros(systemSize^2,numt);


G(occupySite0,1) = lognrnd(G0,sigG0,[1,length(occupySite0)]);
C(occupySite0,1) = lognrnd(C0,sigC0,[1,length(occupySite0)]);
age(occupySite0,1) = datasample([1:1:maxAge],length(occupySite0));%uniform distr.
biom(occupySite0,1)= maxBiom./(1+exp(-k.*(age(occupySite0,1)-satAge)));

beta = 1;   %in growth function 
B    = 0.05;
m0   = 0.001;
eps  = 0.2; %cost of facilitation 

alpha = 0.1;%effect of h on dh/dt
J     = 0.0015;

gamma = 0.1; %what age mortality increases
dd = 0.1;%dispersal probably decline with distance

seedMin = 0.5;%relationship between biomass and seed production
kk = 10;      %relationship between biomass and seed production

p_disp = 0.01; %threshold dispersal probability

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
    pDisp    = 1./exp(dd.*distance);%0.0464
    
    sourceBiom = biom(sourceLoc,t);
    pSeed = seedMin+(1-seedMin)./(1+exp(-kk.*(sourceBiom-satBiom)));
    pSeedM= repmat(pSeed,1,length(sinkLoc));
    
    sourceC    = C(sourceLoc,t);
    sinkH      = h(sinkLoc,t); 
    fitness    = repmat(sourceC,1,length(sinkLoc))-...
        repmat(sinkH',length(sourceLoc),1);
    fitness    = 1./exp(beta.*fitness.^2);% 0.0319
    
    prob = ((pDisp.*pSeedM).*fitness);
    [ii,jj] = max(prob); %jj is row ID of max prob
    newSites = sinkLoc(find(ii>p_disp));
    parentSite = sourceLoc(jj(find(ii>p_disp)));
    
    G(newSites,t+1) = G(parentSite,t);
    C(newSites,t+1) = C(parentSite,t);
    biom(newSites,t)= 0.01;
    occupySite = unique([occupySite;newSites]);
    
    %% growth    
    rGrow = (biom(occupySite,t).*(1-biom(occupySite,t))).*...
        B.*(1 + 1./exp(beta.*(C(occupySite,t)-h(occupySite,t)).^2))...
        -eps.*(biom(occupySite,t).*G(occupySite,t));
    biom(occupySite,t+1) = biom(occupySite,t) + rGrow;
    
    %% landscape change 
    h(:,t+1)=h(:,t)+(J+G(:,t).*biom(:,t))./exp(alpha*h(:,t)./(hmax-h(:,t)));
    
    %% increase the age of the rest of plants 
    age(occupySite,t+1)  = age(occupySite,t)+1;
    
    %% mortality 
    m = m0./exp((maxAge-age(occupySite,t)).*gamma);
    mortLoc = find(m>rand(length(occupySite),1));
    
    biom(occupySite(mortLoc),t+1) = 0;
    G(occupySite(mortLoc),t+1)    = 0;
    C(occupySite(mortLoc),t+1)    = 0;
    age(occupySite(mortLoc),t+1)  = 0;
      
end
toc

histogram(C(find(C(:,1)>0),1));hold on
histogram(C(find(C(:,1000)>0),1000));hold on
histogram(C(find(C(:,2000)>0),2000));hold on

subplot(1,2,1)
h_plot = reshape(h(:,2000),systemSize,systemSize);
surf(meshgrid(1:1:systemSize),meshgrid(1:1:systemSize)',h_plot)
subplot(1,2,2)
biom_plot = reshape(biom(:,2000),systemSize,systemSize);
surf(meshgrid(1:1:systemSize),meshgrid(1:1:systemSize)',biom_plot)


