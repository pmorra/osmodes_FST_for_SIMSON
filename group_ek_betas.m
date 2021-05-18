function [ekbN,bN] = group_ek_betas(Lz,osbeta,ek)
%
% Computes the energy distribution as function of beta only
% NB: it can be used also for alphas, or gammas
%
% INPUT:  Lz: domain length
%         osbeta: betas available in the original distribution
%         ek: energy available in the original distribution
%
% OUTPUT: ekbN: energy per beta
%         bN:   beta
%
% Pierluigi Morra, 2020
%

check_out = false;
nb = 0; c = 1; n = 0;
while ~check_out
  bn = 2*pi/Lz*n;
  n = n+1;
  [check,idfnd] = find(bn == abs(osbeta));
  if isempty(check)
    continue
  else
    ekbN(c) = sum(ek(idfnd));
    bN(c) = bn;
    nb = nb+length(idfnd);
    if nb == length(osbeta)
      check_out = true;
      break
    end
    c = c+1;
  end
end