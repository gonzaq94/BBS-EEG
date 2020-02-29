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
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
title('original source configuration: two source regions','FontSize',18);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%% Gibbs sampler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Shat = Gibbs_sampler(X(:,id),G);
%Shat
%{
[SOut,LambdaOut]= Gibbs_sampler(X(:,id),G);

title('Localization of a signal of interest (epileptic spike) with the Gibbs sampler');
plot(G*SOut);
hold on;
plot(G*S(:,id));
grid on;
legend('G*s (estimated)','G*s (real)', 'Location','South')
print('images/Gibbs/spikeLocalization.png','-dpng')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MNE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
lambda=1;
s_hat_mne = MNE(X,G,lambda);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_hat_mne(:,id));
title(['MNE algorithm with SNR = 1 and Lambda = ' num2str(lambda)],'FontSize',18);
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


% L curve criterion and Discrepancy principle

SNR = 1;

X=Xs+1/sqrt(SNR)*Noise;

lambda_vec = [linspace(1e-10,10,20),linspace(10,1000,80),linspace(1000,1e10,20)];
lambda_vec = logspace(0,4,120);

err_vec = zeros(length(lambda_vec),1);
norm_s = zeros(length(lambda_vec),1);

for i=1:length(lambda_vec)
    
        s_hat_mne = MNE(X,G,lambda_vec(i));
        err_vec(i) = norm(X-G*s_hat_mne);
        norm_s(i) = norm(s_hat_mne);

end
  
% L criterion
val = 156.9; % from the plot
[ d, ix ] = min( abs( norm_s- val) );

figure;
ax = axes;
plot(norm_s,err_vec);
ax.XScale = 'log';
ax.YScale = 'log';
title('L-curve with SNR = 1');
xlabel('Norm of s');
ylabel('MSE(x,Gs)');
grid on;
dim = [.2 .5 .3 .3];
p = [30 100 200 300];
annotation('textbox',dim,'String',['Lambda=' num2str(lambda_vec(ix))],'FitBoxToText','on','Position', p);
print('images/MNE/L-curvev3.png','-dpng');

%Discrepancy principle
%{
noise_pwr = norm(Noise)*ones(length(lambda_vec),1);

[xi,yi] = polyxpoly(lambda_vec,err_vec,lambda_vec,noise_pwr);

figure;
ax = axes;
plot(lambda_vec,err_vec);
hold on;
plot(lambda_vec,noise_pwr);
title('Discrepancy principle with SNR = 1');
xlabel('Lambda');
%xlim([0, 1]);
%ylim([0, norm(Noise)*2]);
ax.XScale = 'log';
ax.YScale = 'log';
grid on;
legend('Reconstruction MSE','Noise power');
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',['Lambda min=' num2str(xi)],'FitBoxToText','on');
%print('images/MNE/Discrepancy.png','-dpng');
%}
% Generalized cross validation
%{
SNR = 1;

X=Xs+1/sqrt(SNR)*Noise;

lambda_vec = [linspace(1e-10,1e-2,10),linspace(1e-2,1,15),linspace(1,400,100),linspace(400,1e10,80)];
GCV = zeros(length(lambda_vec),1);

sizeG = size(G);

for i=1:length(lambda_vec)
    
        s_hat_mne = MNE(X,G,lambda_vec(i));
        GCV(i) = norm(X-G*s_hat_mne)/(trace(eye(size(G,1))-G*transpose(G)*inv(G*transpose(G)+lambda_vec(i)*eye(size(G,1)))))^2;

end
  
[C,i] = min(GCV);
figure;
ax = axes;
plot(lambda_vec,GCV)
title('Generalized cross validation function with SNR = 1')
ylabel('GCV')
xlabel('Lambda')
ax.XScale = 'log';
ax.YScale = 'log';
grid on;
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',['Lambda min=' num2str(lambda_vec(i))],'FitBoxToText','on');
print('images/MNE/GCV2.png','-dpng');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SISSY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
s_hat_sissy = SISSY(X,G,variation_operator(mesh,'face'),1,0.1,60);

figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),s_hat_sissy(:,id));

plot_eeg(s_hat_mne(idx_electrodes,:),max(max(s_hat_mne(idx_electrodes,:))),256,channel_names);
title('cleaned EEG data','FontSize',18);
%}
%immse(S,s_hat_sissy);

% variation of the regularization parameter
%{
lambda_vec = [linspace(0.01,100,6),linspace(100,1000,4) ];
err_vec = zeros(length(lambda_vec),1);

for i=1:length(lambda_vec)
        i
        s_hat_sissy = SISSY(X,G,variation_operator(mesh,'face'),lambda_vec(i),0.1,40);
        err_vec(i) = immse(S,s_hat_sissy);

end
    
plot(lambda_vec,err_vec)
title('Variation of the regularization parameter with SNR = 1')
ylabel('Mean squared error')
xlabel('Lambda')
grid on;
print('images/SISSY/lambdaVariation.png','-dpng');
%}
% variation of alpha
%{
alpha_vec = [linspace(0,0.1,5),linspace(0.1,1,5) ];
err_vec = zeros(length(alpha_vec),1);

for i=1:length(alpha_vec)
        i
        s_hat_sissy = SISSY(X,G,variation_operator(mesh,'face'),200,alpha_vec(i),40);
        err_vec(i) = immse(S,s_hat_sissy);

end
    
plot(alpha_vec,err_vec)
title('Variation of alpha with SNR = 1 and Lambda=200')
ylabel('Mean squared error')
xlabel('Alpha')
grid on;
print('images/SISSY/alphaVariation.png','-dpng');
%}

%%% lambda paramter optimization
%{
%lambda_vec = linspace(1e-2,1e4,5);
lambda_vec = [linspace(0.01,100,6),linspace(100,1000,4) ];
loss_function = zeros(length(lambda_vec),1);

T = variation_operator(mesh,'face');

alpha = 0.1;

for i=1:length(lambda_vec)
    
        i
        s_hat_sissy = SISSY(X,G,T,lambda_vec(i),alpha,40);
        
        threshold = 0.01*max(max(s_hat_sissy));
        
        s_aux = s_hat_sissy;
        s_aux(find(s_aux<threshold))=0;
        
        Ts = T*s_hat_sissy;
        Ts(find(Ts<threshold)) = 0;
        
        loss_function(i) = norm(Ts)+alpha*norm(s_aux);

end
    
plot(lambda_vec,loss_function)
title('Variation of the contraint with lambda for SNR = 1')
xlabel('Lambda');
ylabel('L0(Ts) + alpha L0(s)');
grid on;
print('images/SISSY/lambdaOptimization.png','-dpng');
%}