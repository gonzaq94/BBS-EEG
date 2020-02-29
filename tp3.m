clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source localization algorithms to this time point
[~,id]=max(mean(S,1));

%visualize original source distribution
%figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
%title('original source configuration: two source regions','FontSize',18); axis off;


%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

active_dipoles = find(mean(S,2)>0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Qualitative comparisons %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Different SNRs
%{
%signal to noise ratio
SNR_vec=[0.1 1 10];


for i=1:length(SNR_vec)
    
    %generate noisy data according to given SNR
    X=Xs+1/sqrt(SNR_vec(i))*Noise;
    
    s_hat_mne = MNE(X,G,228);
    s_hat_sissy = SISSY(X,G,variation_operator(mesh,'face'),10,0.1,40);
    
    figure;
    subplot(1,2,1);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_hat_mne(:,id));
    title(['MNE algorithm with SNR = ' num2str(SNR_vec(i))]);
    
    subplot(1,2,2);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_hat_sissy(:,id));
    title(['SISSY algorithm with SNR = ' num2str(SNR_vec(i))]);
    
end
%}
% Correlated noise
%{
% generate random noise, that corresponds to the background activity model
Snoise = zeros(size(S));
Snoise(find(S==0)) = randn(size(find(S==0)));
Noise = G*Snoise;

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');


SNR_vec=[0.1 1 10];


for i=1:length(SNR_vec)
    
    %generate noisy data according to given SNR
    X=Xs+1/sqrt(SNR_vec(i))*Noise;
    
    s_hat_mne = MNE(X,G,228);
    s_hat_sissy = SISSY(X,G,variation_operator(mesh,'face'),10,0.1,40);
    
    figure;
    subplot(1,2,1);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_hat_mne(:,id));
    title(['MNE algorithm with SNR = ' num2str(SNR_vec(i))]);
    
    subplot(1,2,2);
    trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_hat_sissy(:,id));
    title(['SISSY algorithm with SNR = ' num2str(SNR_vec(i))]);
    
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Quantitative comparisons %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
SNR = 1;

% generate dipole positions matrix
i=1;
r_grid = zeros(size(mesh.f));
while(i<=size(mesh.f,1))
    r_grid(i,:)= mean(mesh.v(mesh.f(i,:),:));
    i =i+1;
end
%save('r_grid.mat',r_grid);
inactive_dipoles = find(S==0);
%active_dipoles = find(mean(S,2)>0);
%active_dipoles = find(S(:,id)>4);
active_dipoles = unique(mod(find(abs(S)>3),size(S,1)));
active_dipoles(find(active_dipoles==0)) = size(S,1);

DLE_matrix = zeros(10,2);
threshold = 2;

for i=1:10
    i
    % generate random noise, that corresponds to the background activity model
    Snoise = zeros(size(S));
    Snoise(inactive_dipoles) = randn(size(inactive_dipoles));
    Noise = G*Snoise;
    Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');
    
    X=Xs+1/sqrt(SNR)*Noise;
    
    s_hat_mne = MNE(X,G,228);
    s_hat_sissy = SISSY(X,G,variation_operator(mesh,'face'),10,0.1,40);
    
    active_est_dipoles_sissy = unique(mod(find(abs(s_hat_sissy)>threshold),size(s_hat_sissy,1)));
    active_est_dipoles_sissy(find(active_est_dipoles_sissy==0)) = size(s_hat_sissy,1);
    
    active_est_dipoles_mne = unique(mod(find(abs(s_hat_mne)>threshold),size(s_hat_mne,1)));
    active_est_dipoles_mne(find(active_est_dipoles_mne==0)) = size(s_hat_mne,1);

    DLE_matrix(i,1) = DLE(active_dipoles,active_est_dipoles_mne,r_grid);
    DLE_matrix(i,2) = DLE(active_dipoles,active_est_dipoles_sissy,r_grid);
    
end

boxplot(DLE_matrix,'Labels',{'MNE','SISSY'});
grid on;
title('DLE boxplot for SNR = 1');
grid on;
print('images/Performance Analysis/DLEboxplot.png','-dpng');

c1=zeros(size(S(:,1)));
c1(active_dipoles)=1;
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),c1);

c2=zeros(size(S(:,1)));
c2(active_est_dipoles_mne)=1;
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),c2);4
%}

%%% DLE vs SNR plot

SNR_vec = linspace(0.1,10,5);

% generate dipole positions matrix
i=1;
r_grid = zeros(size(mesh.f));
while(i<=size(mesh.f,1))
    r_grid(i,:)= mean(mesh.v(mesh.f(i,:),:));
    i =i+1;
end
%save('r_grid.mat',r_grid);
inactive_dipoles = find(S==0);
%active_dipoles = find(mean(S,2)>0);
%active_dipoles = find(S(:,id)>4);
active_dipoles = unique(mod(find(abs(S)>3),size(S,1)));
active_dipoles(find(active_dipoles==0)) = size(S,1);

DLE_matrix = zeros(length(SNR_vec),2);
threshold = 2;
Nit=10;

for i=1:Nit
    i
    % generate random noise, that corresponds to the background activity model
    Snoise = zeros(size(S));
    Snoise(inactive_dipoles) = randn(size(inactive_dipoles));
    Noise = G*Snoise;
    Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');
    
    for j=1:length(SNR_vec)
    
        X=Xs+1/sqrt(SNR_vec(j))*Noise;

        s_hat_mne = MNE(X,G,228);
        s_hat_sissy = SISSY(X,G,variation_operator(mesh,'face'),10,0.1,40);

        active_est_dipoles_sissy = unique(mod(find(abs(s_hat_sissy)>threshold),size(s_hat_sissy,1)));
        active_est_dipoles_sissy(find(active_est_dipoles_sissy==0)) = size(s_hat_sissy,1);

        active_est_dipoles_mne = unique(mod(find(abs(s_hat_mne)>threshold),size(s_hat_mne,1)));
        active_est_dipoles_mne(find(active_est_dipoles_mne==0)) = size(s_hat_mne,1);

        DLE_matrix(j,1) = DLE_matrix(j,1) + DLE(active_dipoles,active_est_dipoles_mne,r_grid);
        DLE_matrix(j,2) = DLE_matrix(j,2) + DLE(active_dipoles,active_est_dipoles_sissy,r_grid);
        
    end
       
end

DLE_matrix = DLE_matrix /Nit;

plot(SNR_vec,DLE_matrix(:,1));
hold on;
plot(SNR_vec,DLE_matrix(:,2));
grid on;
title('DLE boxplot for SNR = 1');
xlabel('SNR');
ylabel('Mean DLE');
title('Mean DLE for varying SNR')
legend('MNE','SISSY');
print('images/Performance Analysis/DLEvsSNR.png','-dpng');