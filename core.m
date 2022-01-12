% This code is the topic model based on the von Mises-Fisher (vMF) distribution
% See paper for more details
% All experiments were done using MATLAB R2013b

function [ clust,muMat ] = core(n_UnitVectors, mNm, nm, m_apprKappa, K, params)
rng('default');
addpath('utils/');
EM_MAX_ITER = params.EM_MAX_ITER;
EM_CONVERGE = params.EM_CONVERGE;
EM_PERCENT_ELEM = params.EM_PERCENT_ELEM;
saveVars = params.saveVars;
nameVars = params.nameVars;

[N,dim] = size(n_UnitVectors);
[M,~] = size(mNm);

%% Create plate notation
M = length(unique(nm(:,2)));
plate = {};
for m=1:M
    mb = find(nm(:,2)==m);
    plate{m} = nm(mb,1);
end;

% Initialization
[ muMat,alpha,gamma_mk,phiMat ] = initVals(n_UnitVectors,M,N,K,params.methodmu);

mainIter = 1;
toContinue = 1;
converged = 0;
completell_old = 100;
T = 1000;
alpha_suffstats = 0;

n_mKappa = m_apprKappa(nm(:,2),2);
kappa = repmat(n_mKappa,1,K);
while(toContinue == 1)

    gamma_mk_old = gamma_mk;
    for m=1:M
        sumphiMat(m,:) = sum(phiMat(plate{m},:),1);
    end;
    gamma_mk = alpha + sumphiMat;

    phiMat_old = phiMat;
    psigamma = psi(gamma_mk);
    psisumgamma = psi(sum(gamma_mk,2));
		  
	logNormalize  = (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(dim/2-1,kappa);
	logProbMat    = ( bsxfun(@times,(n_UnitVectors * muMat'),kappa) + logNormalize )/T;
	  
    [~,clust] = max(logProbMat,[],2);
    
    for m=1:M	
    	Nm = mNm(m,2);
    	for n=1:Nm
    		onumb = plate{m}(n);	
    		for h=1:K	
				phiMat(onumb,h) = exp(psigamma(m,h) - psisumgamma(m) + logProbMat(onumb,h));
    		end;
            phiMatCopy(onumb,:) = phiMat(onumb,:);
			phiMat(onumb,:) = bsxfun(@rdivide,phiMat(onumb,:),sum(phiMat(onumb,:)));
    	end;
        
		gamma_sum = sum(gamma_mk(m,:),2);
		alpha_suffstats = sum(psi(gamma_mk(m,:)),2);
		alpha_suffstats = alpha_suffstats - K * psi(gamma_sum);    
    end;
    
    SumphiMat = sum(phiMatCopy,1);
    mixprop = bsxfun(@rdivide,SumphiMat,sum(SumphiMat));	
	mu_old = muMat;
	
	diff      = 1;
	epsilon   = 0.0001;
	value     = 100;
	iteration = 1;
	while (diff > epsilon)
		iteration = iteration + 1;
		oldvalue  = value;

		logNormalize2  = repmat(log(mixprop),N,1) + (dim/2-1)*log(kappa) - (dim/2)*log(2*pi) - logbesseli(dim/2-1,kappa);
		logProbMat2    = ( bsxfun(@times,(n_UnitVectors * muMat'),kappa) + logNormalize2 )/T;	

		logSum2        = log(sum(exp(logProbMat2),2));					       
		logProbMat2    = logProbMat2 - logSum2*ones(1,K);  
		value = sum(sum(logProbMat2));
		muMat     = exp(logProbMat2')*n_UnitVectors;
  
		for h=1:K
			normMu   = sqrt(muMat(h,:)*muMat(h,:)');   
			muMat(h,:)  = muMat(h,:)/normMu;
		end;
		diff = abs(value - oldvalue);
	end;
	
	optalpha = opt_alpha(alpha_suffstats, M, K);
	alpha = optalpha;
	
	diff_gamma_mk = abs(gamma_mk_old - gamma_mk);
	diff_phiMat = abs(phiMat_old - phiMat);
	diff_mu = abs(mu_old - muMat);
    
    cnt_diffgamma = length(find(diff_gamma_mk < EM_CONVERGE));
	[row,col] = size(gamma_mk);
	percentDiffgamma = cnt_diffgamma/(row * col) * 100;
	
    cnt_diffphi = length(find(diff_phiMat < EM_CONVERGE));
	[row,col] = size(phiMat);
	percentDiffphi = cnt_diffphi/(row * col) * 100;
	
    cnt_diffMu = length(find(diff_mu < EM_CONVERGE));
	[row,col] = size(muMat);
	percentDiffMu = cnt_diffMu/(row * col) * 100;
	
    statPass = 0;
    if ((percentDiffphi >= EM_PERCENT_ELEM &...
       percentDiffMu >= EM_PERCENT_ELEM) | (mainIter == EM_MAX_ITER))
       statPass = 1;
    end;
    
	if (statPass == 1)
        if (length(unique(clust)) == K)
            toContinue = 0;
            break;
        end;
        mainIter = mainIter + 1;
    else
        mainIter = mainIter + 1;
	end;
    
    if (mainIter == EM_MAX_ITER)
        break;
    end;    
end;

if (saveVars==1)
    save(nameVars);
end;


function optalpha = opt_alpha(alpha_suffstats, M, K)
	a = 0; log_a = 0; init_a = 100;
	f = 0; df = 0; d2f = 0;
	iter = 0;
	log_a = log(init_a);
	
	NEWTON_THRESH = 1e-5;
	MAX_ALPHA_ITER = 1000;
	
    toContinue = true;
	while (toContinue)
		iter = iter + 1;
		a = exp(log_a);
		if (isnan(a))
            init_a = init_a * 10;            
            a = init_a;
            log_a = log(a);
        end;
        df = d_alhood(a, alpha_suffstats, M, K);
        d2f = d2_alhood(a, M, K);
        log_a = log_a - df/(d2f * a + df);	
        
        if ((abs(df)> NEWTON_THRESH) & (iter < MAX_ALPHA_ITER))
            toContinue = true;
        else
            toContinue = false;
        end;
	end;
	optalpha = exp(log_a);	
return

function df = d_alhood(a, ss, M, K)
    df = M * (psi(K * a) - K * psi(a)) + M * ss;
return

function d2f = d2_alhood(a, M, K)
	d2f = M * (K * K * psi(1,(K * a)) - K * psi(1,(a)));
return

