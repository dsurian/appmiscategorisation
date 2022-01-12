function [ muMat,alpha,gamma_mk,phiMat ] = initVals(n_UnitVectors,M,N,K,methodmu)

muMat = initMu(n_UnitVectors,K,methodmu);

initVal_alpha = 1;
initVal_gamma = 1;
initVal_phi = 1;
alpha = initVal_alpha;
gamma_mk = ones(M,K) * initVal_gamma;
phiMat = ones(N,K) * initVal_phi;
