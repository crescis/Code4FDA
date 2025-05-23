function fdobj_reco = reco(pca,pc,Ntrunc)
% Reconstruction of multivariate hydrographic profiles
%
% RECO This function reconstructs multivariate hydrographic profiles with 
% a chosen number of Principal Components (PCs).
%
% ARGUMENTS
% PCA ... structure containing the vertical modes to project on
% PC  ... Principal Components of the multivariate profiles, computed with
% the function reco.
% NTRUNC number of PCs to use in the reconstruction, default is set to the total number of PC, \code{te = mybn}.
% @return \code{recotemp} and \code{recosali} matrix containing the reconstruction of \code{temp} and \code{sali} with the number of PCs \code{te}.
%
% RETURN
% FDOBJ_RECO ... An FD object containing the reconstructed multivariate
% profiles. They can be evaluated at selected levels Pi with the function eval_fd.
%
% DEPENDENCIES
% The method uses the fdaM Toolbox by Jim Ramsay.
% http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/
% You will need to install this toolbox and add it to the matlab path to use this software
%
% CONTACT
% This code was written by Etienne Pauthenet, David Nerini and Fabien Roquet. 
% Questions, comments and bugs can be sent to: 
% etienne.pauthenet@gmail.com
% 
% REFERENCES 
% Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1
% Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
%

basis = pca.basis;
prange = getbasisrange(basis);
nbas = pca.nbas;
ndim = pca.ndim;

if ~exist('Ntrunc','var'), Ntrunc = ndim*nbas; end

nobs = size(pc,1);
coef = zeros(nbas,nobs,ndim);
for kk=1:ndim,
    coef(:,:,kk) = repmat(pca.Cm((kk-1)*nbas+1:kk*nbas)',1,nobs) + pca.vectors((kk-1)*nbas+1:kk*nbas,1:Ntrunc)*pc(:,1:Ntrunc)';
end
%T_domain = (0:(50-1))/(50-1); 
%fdobj_reco = data2fd(T_domain',coef,basis,2,0.0001);
fdobj_reco = fd(coef,basis,pca.fdnames);
