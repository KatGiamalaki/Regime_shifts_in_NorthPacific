% Computation of theta and d for each days ih the database. 
% The dataset should be a .mat file of size (number_of_time_observations*lon*lat)
% Code provided by Davide Faranda based on Faranda et al.,(2017) Dynamical proxies of North Atlantic predictability and extremes. Sci. Rep. 7, 41278, doi:10.1038/srep41278

% % Load data from Data_Prep.m

load('NP_slp');

aa = NP_slp;  % current shape is lon*lat*n_days
aa = permute(aa,[3 1 2]); 
aa = reshape(aa, [size(aa,1) size(aa,2)*size(aa,3)]); % reshape to n_days lon*lat

clearvars -except aa

%%
% Dynamical Systems Analysis

tic
quanti = 0.98;

a = 1;
for j = 1:size(aa,1)
    j   
        distance = pdist2(aa(j,:),aa);
        logdista = -log(distance);
        theta(j) = extremal_Sueveges(logdista,quanti);
        thresh = quantile(logdista, quanti);
        logdista = sort(logdista);
        findidx = find(logdista > thresh,1);
        logextr = logdista(findidx:end-1);

        [tpar, tpari] = gpfit(logextr - thresh);
        Acsi(j) = tpar(1);
        Asigma(j) = tpar(2);
        icsi(j) = tpari(1,1) - Acsi(j); 
        a = a + 1;    
end 

toc


