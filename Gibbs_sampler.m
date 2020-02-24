function [SOut,LambdaOut]=Gibbs_sampler(X,A)

%number of iterations
Niter=100;%100
epsilon = 0.001;

%get number of sensors and number of dipoles
[N,D]=size(A);

% constants 
nA = sum(A.^2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% implement the Gibbs sampler here according to the pseudocode
%
% store vectors q and s of each iteration in matrices Q (D x Niter) and S
% (D x Niter)
% use variables sigma_n2 and sigma_s2 for the variances of noise and
% signals

Q = zeros(D,Niter);
S = zeros(D,Niter);

% initialization
q = Q(:,1);
s = S(:,1);
lambda = betarnd(1,1);
%sigma_n2 = 1/gamrnd(1,1e4);
%sigma_s2 = 1/gamrnd(1,1);

sigma_n2 = InverseGamma(1,1e-4);
sigma_s2 = InverseGamma(1,1);

for it=1:Niter
    %step 1: sample (qi,si)
    e = X-A*s;
    for i=1:N

       ei = e+A(:,i)*s(i); %A(:,i) is the i-th column of A
       sigma_i = sigma_n2*sigma_s2/(sigma_n2+sigma_s2*norm(A(:,i))^2);
       mu_i = sigma_i/sigma_n2*A(:,1)'*ei;
       v_i = lambda*(sqrt(sigma_i)/sqrt(sigma_s2))*exp(mu_i^2/(2*sigma_i)); 
       lambda_i = v_i/(v_i+1-lambda);

       if isnan(lambda_i)
          
           qi = 1;
       else
           qi = binornd(1,lambda_i);
       end

       if qi == 0
           si = 0;
       else
           si = mu_i+sqrt(sigma_i)*randn();
       end
       
       e = ei - A(:,i)*si; %update e

       Q(i,it) = qi;
       S(i,it) = si;
    end
    
    q = Q(:,it);
    s = S(:,it);
    %step 2: sample sigma_n2
    %sigma_n2 = 1/gamrnd(N/2,2./norm(X-A*s)^2); %things missing
    sigma_n2 = InverseGamma(N/2+epsilon, norm(X-A*s)^2/2+epsilon);
    %step 3: sample lambda
    L = sum(q);
    lambda = betarnd(1+L,1+D-L);
    %step 4: sample sigma_s2
    %sigma_s2 = 1/gamrnd(L/2,2./norm(s)^2); %things missing
    sigma_s2 = InverseGamma(L/2+epsilon, norm(s)^2/2+epsilon);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Make the estimation from Q and S using the MAP criterium
%Q(:,Niter/2:Niter)

for j = 1:Niter/2
    q = Q(:,Niter/2+j); 
    idx = find(q==1); 
    R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2; 
    S_opt = (R\(A(:,idx)'*X))/sigma_n2; 
    cout(j) = - S_opt'*R*S_opt/sigma_n2; 
end

[Val, OptIdx] = min(cout); 
q = Q(:,Niter/2+OptIdx);
idx = find(q==1); 
R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2; 
SOut = zeros(D,1); 
SOut(idx) = (R\(A(:,idx)'*X))/sigma_n2;
LambdaOut = length(idx)/D; 



