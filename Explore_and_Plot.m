%% Load data resulting from Dynamical_Proxies.m and Data_Prep.m

load('NP_slp.mat')
load('Dyn_prox.mat')
load('time_full.mat') % load time array in format dd, yyy, mm
time_obs = time_full'; 

%% Calculate delta(z) and quantiles for d and th

deltaz = 1./Asigma;
Y = deltaz;
p = 0.98;
Z98 = quantile(Y,p);
p = 0.02;
Z02 = quantile(Y,p);

X = theta;
r = 0.02;
T02 = quantile(X,r);
r = 0.98;
T98 = quantile(X,r);

D = nanmean(deltaz);
T = nanmean(theta);

%% Make vectors
theta = theta(:,1:size(aa,1));
th1 = theta(1:size(aa,1));
th1(2,:) = time_obs(2,:);
th1(3,:) = 1:size(theta,2);
th1(4,:) = time_obs(3,:);

deltaz = deltaz(:,1:size(aa,1));

del1 = deltaz(1:size(aa,1));
del1(2,:) = time_obs(2,:);
del1(3,:) = 1:size(theta,1);
del1(4,:) = time_obs(3,:);

%% Isolating 1976-1977 data for major regime shift

II = th1(2,:)==1977;
th77 = th1(:,II);
II = th1(2,:)==1976;
th76 = th1(:,II);

%% Defining theta and delta thresholds  

clearvars TH TL ZH ZL Tmax Tmin Zmax Zmin

        TH = find(th1(1,:)>=T98 & del1(1,:)<=Z98 & del1(1,:)>=Z02);
        Tmax(1,:) = th1(1,TH);
        Tmax(2,:) = th1(2,TH);
        Tmax(3,:) = deltaz(1,TH);
        Tmax(4,:) = time_obs(1,TH);
        Tmax(5,:) = time_obs(3,TH);
        
SLPTmax = NP_slp(:,:,TH);

        
        TL = find(th1(1,:)<=T02 & del1(1,:)<=Z98 & del1(1,:)>=Z02);
        Tmin(1,:) = th1(1,TL);
        Tmin(2,:) = th1(2,TL);
        Tmin(3,:) = deltaz(1,TL);
        Tmin(4,:) = time_obs(1,TL);
        Tmin(5,:) = time_obs(3,TL);
        
SLPTmin = NP_slp(:,:,TL);
 
        ZH = find(del1(1,:)>=Z98 & th1(1,:)>=T02 & th1(1,:)<=T98);
        Zmax(1,:) = th1(1,ZH);
        Zmax(2,:) = th1(2,ZH);
        Zmax(3,:) = deltaz(1,ZH);
        Zmax(4,:) = time_obs(1,ZH);
        Zmax(5,:) = time_obs(3,ZH);
        
SLPZmax = NP_slp(:,:,ZH);
    
        ZL = find(del1(1,:)<=Z02 & th1(1,:)>=T02 & th1(1,:)<=T98);
        Zmin(1,:) = th1(1,ZL);
        Zmin(2,:) = th1(2,ZL);
        Zmin(3,:) = deltaz(1,ZL);
        Zmin(4,:) = time_obs(1,ZL);
        Zmin(5,:) = time_obs(3,ZL);
        
SLPZmin = NP_slp(:,:,ZL);

Zmin2 = Zmin;
%% Single to double 

SLPTmax = double(SLPTmax);
SLPTmin = double(SLPTmin);
SLPZmax = double(SLPZmax);
SLPZmin = double(SLPZmin);

%% Plot delta z against theta ++ years of interest highlighted 

x1 = Z02;
x2 = Z98;
y3 = T02;
y4 = T98;

figure;
subplot(3,3,[4 5 7 8]);
scatter(deltaz,theta,'.','m'); xlabel('d ( \zeta )'); ylabel('\theta ( \zeta )'); 
y1 = get(gca,'ylim'); x3=get(gca,'xlim');
ylabel('\theta ( \zeta )','FontSize', 12, 'FontWeight', 'bold') ; 
xlabel('d(\zeta)','FontSize', 12, 'FontWeight', 'bold');
box on
hold on
plot([x1 x1],y1,'LineStyle','--','LineWidth',2, 'Color',[0 0 0]);
hold on
plot([x2 x2],y1,'LineStyle','--','LineWidth',2, 'Color',[0 0 0]);
hold on 
plot(x3, [y3 y3],'LineStyle','--','LineWidth',2, 'Color',[0 0 0]);
hold on 
plot(x3,[y4 y4],'LineStyle','--','LineWidth',2, 'Color',[0 0 0]);
hold off

subplot(3,3,2);
pcolor(lon,lat,squeeze(nanmean(SLPTmax,3)));shading interp;xlabel('Longitude');ylabel('Latitude'); title('\theta(\zeta) maximum');fillmap;caxis([99027.77108433735,103968.2469879518]);%title(colorbar,'Pa');
subplot(3,3,9);
pcolor(lon,lat,squeeze(nanmean(SLPTmin,3)));shading interp;xlabel('Longitude');ylabel('Latitude');title('\theta(\zeta) minimum');fillmap;caxis([99027.77108433735,103968.2469879518]);%caxis([-0.5 0.5]);%title(colorbar,'Pa');
subplot(3,3,6);
pcolor(lon,lat,squeeze(nanmean(SLPZmax,3)));shading interp; xlabel('Longitude');ylabel('Latitude');title('d(\zeta) maximum');fillmap;caxis([99027.77108433735,103968.2469879518]);%caxis([-0.5 0.5]);%colorbar;title(colorbar,'Pa');
subplot(3,3,1);
pcolor(lon,lat,squeeze(nanmean(SLPZmin,3)));shading interp;xlabel('Longitude');ylabel('Latitude');title('d(\zeta) minimum');fillmap;caxis([99027.77108433735,103968.2469879518]);%caxis([-0.5 0.5]);title('\delta(\zeta) minimum');%colorbar;title(colorbar,'Pa');

%%
clearvars nelements xcentres alpha a b ii II count_year count_month average_days average_month

alpha = Zmin;
a = 1948;
b = 1;
for ii = 1:size(alpha,2)
    II = find(alpha(2,:)==a);
    count_year(b,1) = length(II);
    a = a+1;
    b = b+1;
end

count_year = count_year(1:68,1)';
average_days = nanmean(count_year);

a = 1;
b = 1;
for ii = 1:size(alpha,2)
    II = find(alpha(5,:)==a);
    count_month(b,1) = length(II);
    a = a+1;
    b = b+1;
end
count_month = count_month(1:12,:)';
average_monthd = nanmean(count_month);
%%
clearvars nelements xcentres alpha a b ii II count_year count_month average_days average_month

alpha = Zmin;
a = 1948;
b = 1;
for ii = 1:size(alpha,2)
    II = find(alpha(2,:)==a);
    count_year(b,1) = length(II);
    a = a+1;
    b = b+1;
end

count_year = count_year(1:68,1)';
average_days = nanmean(count_year);

a = 1;
b = 1;
for ii = 1:size(alpha,2)
    II = find(alpha(5,:)==a);
    count_month(b,1) = length(II);
    a = a+1;
    b = b+1;
end
count_month = count_month(1:12,:)';
average_monthd = nanmean(count_month);

% calculate the high quantile to select the most persistent years?
% use counts per year

q_count = quantile(count_year(1,:),0.98);
EE = find(count_year(1,:)>=q_count);
y3 = q_count;
xcentresy = 1:68;
nelementsm = count_year(1,:);

figure;
hb = bar(xcentresy,nelementsm); xlabel('Years','FontSize', 12, 'FontWeight', 'bold');
ylabel('No of occurrences','FontSize', 12, 'FontWeight', 'bold'); 
xlim([0 68]);
set(hb,'FaceColor',[.2 .2 .2]);
y1 = get(gca,'ylim'); x3 = get(gca,'xlim');
hold on
plot(x3, [y3 y3],'LineStyle','--','LineWidth',2, 'Color',[0 0 0]);
hold off

%%
% Calculate the high quantile to select the most persistent years
% Use theta(z)

q_th = quantile(Zmin2(1,:),0.98); 
EE_th = find(Zmin2(1,:)>=q_th);

y3 = q_th;
figure; plot(Zmin2(2,:),Zmin2(1,:),'k.');
pro1 = Zmin2(2,EE_th)';
labels = char(num2str(pro1));
text(Zmin2(2,EE_th),Zmin2(1,EE_th),labels);
y1 = get(gca,'ylim');x3=get(gca,'xlim');
hold on
plot(x3, [y3 y3],'LineStyle','--','LineWidth',2, 'Color',[0 0 0]);
hold off

% clearvars -except del1 lat lon NPslp016 SLPZmin th1 time_sst Zmin