function [ek,kk,nmodes,osalfa,osbeta,osgamma,ymax] = ek_osmodes(fname,etastart)
%
% Computes E(k)*dk spectrum from the osmodes file for SIMSON
%
% NB:  E(k)*dk = || (u,v,w) ||^2  (energy norm)
%
% NB:  E(k)*dk mismatch with theoretical E(k) because the number
%      of modes per shell should be included (max is 10 with a
%      dodecahedron)
%
% INPUT:  fname:    osmodes file name
%         etastart: start of blending function (used for free-stream BC)
%         
%
% OUTPUT: ek:     E(k)*dk
%         kk:     kappas
%         nmodes: number of osmodes
%         osalfa: alphas
%         osbeta: betas
%         osgamma: gammas
%         ymax:   Ly of the domain
%
% Pierluigi Morra, 2020
%

[~,eigf2,eigf3,eigf4,osalfa,osbeta,osgamma,osomega,ymax] = read_osmodes(fname);
nmodes = length(osbeta);
ny = size(eigf2,1);

y = clencurt(ny-1); y = (y+1)*ymax/2;
ek = zeros(nmodes,1);
kk = ek; 

for ii = 1:length(osbeta)
  ek(ii) = en_kin(etastart,y,real(osalfa(ii)),osbeta(ii),...
                      eigf2(:,ii),eigf3(:,ii),eigf4(:,ii));
  k2 = real(osalfa(ii)).^2 + real(osbeta(ii)).^2 + real(osgamma(ii)).^2;
  kk(ii) = sqrt(k2);
end
