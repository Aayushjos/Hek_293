%%%%%%%%%%% SVD / EIG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recon_sharpness = funcAutoFocusSVD(I,kappa)
[M,N] = size(I);
I_norm = I./sqrt(sum(I(:).^2));
I_mean = mean(I_norm(:));
G = I_norm - I_mean;
G_g = G*G';
% [~,Sigma,~] = svd(G_g);
% eigenMat = Sigma'*Sigma;
% diagEigen = diag(eigenMat);
diagEigen = eig(G_g);

recon_sharpness = sum(diagEigen(1:end-kappa));
