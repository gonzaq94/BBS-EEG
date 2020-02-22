clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point
[~,id]=max(mean(S,1));

%visualize original source distribution
%figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
%title('original source configuration: two source regions','FontSize',18); axis off;


%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

%signal to noise ratio
SNR=1;

%generate noisy data according to given SNR
X=Xs+1/sqrt(SNR)*Noise;

%visualize data (for a reduced number of sensors whose indices are 
%specified by idx_electrodes)
%plot_eeg(X(idx_electrodes,:),max(max(X(idx_electrodes,:))),256,channel_names);
%title('noisy EEG data','FontSize',18);

%% Gibbs sampler

%Shat = Gibbs_sampler(X(:,id),A);

%% MNE
%{
s_hat_mne = MNE(X,G,1);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_hat_mne(:,id));

plot_eeg(s_hat_mne(idx_electrodes,:),max(max(s_hat_mne(idx_electrodes,:))),256,channel_names);
title('cleaned EEG data','FontSize',18);

err_mne = immse(S,s_hat_mne)
%}
% variation of the regularization parameter
%{
lambda_vec = linspace(0.1,1000,100);
err_vec = zeros(length(lambda_vec),1);

for i=1:length(lambda_vec)
    
        s_hat_mne = MNE(X,G,lambda_vec(i));
        err_vec(i) = immse(S,s_hat_mne);

end
    
plot(lambda_vec,err_vec)
title('Variation of the regularization parameter with SNR = 1')
ylabel('Mean squared error')
xlabel('Lambda')
grid on;
%}

% variation of the SNR and the regularization parameter
%{
snr_vec = linspace(0.1,10,10);
lambda_vec = linspace(0.1,100,40);
err_vec = zeros(length(lambda_vec),length(snr_vec));

for i=1:length(snr_vec)
    
    for j=1:length(lambda_vec)
        
        X=Xs+1/sqrt(snr_vec(i))*Noise;
        s_hat_mne = MNE(X,G,lambda_vec(j));
        err_vec(j,i) = immse(S,s_hat_mne);
        
    end
    
    plot(lambda_vec,err_vec(:,i), 'DisplayName', ['SNR=' num2str(snr_vec(i))]);
    hold on;

    
end


title('Variation of the regularization parameter with varying SNR')
ylabel('Mean squared error')
xlabel('Lambda')
grid on;
legend('show');
print('images/lambdaVsSNR1.png','-dpng');
%}


% L curve criterion
%{
SNR = 1;

X=Xs+1/sqrt(SNR)*Noise;

lambda_vec = logspace(0.1,10,100);
err_vec = zeros(length(lambda_vec),1);

for i=1:length(lambda_vec)
    
        s_hat_mne = MNE(X,G,lambda_vec(i));
        err_vec(i) = immse(S,s_hat_mne);

end
    
loglog(lambda_vec,err_vec)
title('L-curve with SNR = 1')
ylabel('Mean squared error')
xlabel('Lambda')
grid on;
print('images/L-curve.png','-dpng');

%}

% Generalized cross validation

SNR = 1;

X=Xs+1/sqrt(SNR)*Noise;

lambda_vec = linspace(0.1,10,60);
GCV = zeros(length(lambda_vec),1);

sizeG = size(G);

for i=1:length(lambda_vec)
    
        s_hat_mne = MNE(X,G,lambda_vec(i));
        GCV(i) = immse(Xs,G*s_hat_mne)/(trace(eye()-G*transpose(G)*(G*transpose(G)+lambda_vec(i)*eye(sizeG(1)))))^2;

end
    
plot(lambda_vec,GCV)
title('Generalized cross validation function with SNR = 1')
ylabel('GCV')
xlabel('Lambda')
grid on;
print('images/GCV.png','-dpng');
