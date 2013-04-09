function [results,chain,s2chain]=dramrun(model,data,params,options)
%DRAMRUN  Metropolis-Hastings MCMC run with adaptive delayed rejection (DRAM)
%
% This function generates MCMC chain using DRAM adaptation for a model defined
% by user supplied sum-of-squares function and with additive i.i.d. Gaussian
% errors for the observations. The error variance sigma2 is updated
% using conjugate inverse gamma distribution.
%
% [results,chain,s2chain]=dramrun(model,data,params,options)
%
% input:
%
% model.ssfun =    ; % sum-of-squares function, ss=ssfun(par,data),
%                    % that returns  -2*log(p(y|par))
% model.priorfun = ; % prior "sum-of-squares", priorfun(par,params),
%                    % that returns -2*log(p(par)),
%                    % default: inline('0','x','params')
%
% data   = ;         % extra argument for ssfun (to pass the data etc.)
%
% params.par0   =  ; % initial parameter vector (a row vector)
% params.sigma2 =  1;% initial/prior value for the Gaussian error variance
% params.n0     = -1;% precision of sigma2 as imaginative observations
%                    %   if n0<0, no sigma2 update
% params.n      = ;  % number of actual observations (for sigma2 update)
% params.bounds = ;  % 2*npar matrix of parameter bounds
%                    % default: [-Inf,Inf]
%
% options.nsimu  = 2000;   % length of the chain
% options.qcov   = ;       % proposal covariance matrix
%
% parameters for DRAM
% options.adaptint = 10;  % how often to adapt, if zero, no adaptation
% options.drscale  = 3;   % scale for the second proposal, if zero, no DR
%
% output:
%
% results  structure that contains some info about the run
% chain    nsimu*npar MCMC chain
% s2chain  sigma² chain (if generated)


% calls covupd.m for covariance update and (optionally) gammar_mt.m for
% gamma variates

% this is a 'simple' version for demonstration and educational purposes

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.0 $  $Date: $

% Modified by me to add ss score to the chain, formatting
% and added the multiple DR version


  %% get values from the input structs
  nsimu  = getpar(options,'nsimu',10000);
  % initial parameter vector
  par0   = getpar(params,'par0');
  par0   = par0(:)'; % row vector
  % number of parameters
  npar   = length(par0);
  % 2*npar matrix of parameter bounds
  bounds = getpar(params,'bounds',(ones(npar,2)*diag([-Inf,Inf]))');
  % sum-of-squares function, ssfun(par,data),  -2*log(p(y|theta))
  ssfun  = getpar(model,'ssfun');
  % prior "sum-of-squares", -2*log(p(theta))
  %priorfun = getpar(model,'priorfun',inline('0','x','params'));
  if (isfield(model, 'priorfun'))
    warning('priors not handle in the multiple DR');
  end

  %%% parameters for DRAM
  % how often to adapt, if zero, no adaptation
  adaptint = getpar(options,'adaptint',100);
  % scale for the second proposal, if zero, no DR
  drscale  = getpar(options,'drscale',3);
  % scale for adapting the propsal
  adascale = getpar(options,'adascale',2.4/sqrt(npar));
  % number of trials
  ndelays  = getpar(options,'ndelays',2);
  % log file
  log_file = getpar(options,'log_file', '');

  stall_thresh = getpar(options, 'stall_thresh', 0.1);

  % precision of sigma2 as imaginative observations
  %  if n0<0, no sigma2 update
  %n0  = getpar(params,'n0',-1);
  % number of observations (needed for sigma2 update)
  %if (n0>=0)
  %n = getpar(params,'n');
  %end

  qcov = getpar(options,'qcov'); % proposal covariance

  % blow factor for covariace update
  qcovadj  = getpar(options,'qcovadj',1e-5*min(qcov(qcov~=0)));

  % to DR or not to DR
  dodr = ~(drscale<=0);
  if (~dodr)
    ndelays = 1;
  end

  printint  = getpar(options,'printint',500);
  do_log = ~(isempty(log_file));

  % Multi delays setup
  Rs = cell(ndelays, 1);
  iRs = cell(ndelays, 1);
  Ns = cell(ndelays);
  Ds = cell(ndelays);
  scores = NaN(ndelays+1, 1);
  pars = NaN(ndelays+1, npar);

  R       = chol(qcov); % *adascale; % Cholesky factor of proposal covariance

  for i=1:ndelays
    Rs{i} = R./(drscale^(i-1)); % second proposal for DR try
    iRs{i} = inv(Rs{i});
  end
  chain   = zeros(nsimu,npar+1);  % we store the chain here

  %%%% Assuming gaussian distributions and integrating variance out !!!
  %s20 = 0;
  %if (n0>=0)
  %  s2chain = zeros(nsimu,1);   % the sigma2 chain
  %  s20 = sigma2;
  %else
    s2chain = [];
  %end

  oldpar       = par0(:)';                % first row of the chain
  oldss        = feval(ssfun,oldpar,data);% first sum-of-squares
  %oldprior     = feval(priorfun,oldpar,params);
  acce         = 1;                       %  how many accepted moves
  chain(1,:)   = [oldpar oldss];

%  if (s20>0)
%    s2chain(1,:) = sigma2;
%  end

  % Unique identifier for multiple writers in the same file
  uuid = ['DRAM' num2str(round(rand(1)*100)) ' '];

  if (do_log)
    [fid, err] = fopen(['.' filesep log_file '.dat'], 'a');
    fprintf(fid, [uuid '%% columns="iteration, evalutation | lastbest" (' num2str(clock, '%d/%02d/%d %d:%d:%2.2f') ')\n']);
    fprintf(fid, [uuid '1 : %e |'], oldss);
    fprintf(fid, ' %f',oldpar);
    fprintf(fid, '\n');
    fclose(fid);
  end

  % covariance update uses these to store previous values
  chaincov = [];
  chainmean = [];
  wsum = [];
  lasti = 0;
  if (stall_thresh < 1)
    stall_thresh = nsimu*stall_thresh;
  end
  newi = 0;
  ntot = 0;

  %%% the simulation loop
  for isimu=2:nsimu
%  keyboard


    if (isimu/printint == fix(isimu/printint)) % info on every printint iteration
      fprintf('isimu=%d, %d%% done, accepted: %d%%(%d%%)\n',...
              isimu,fix(isimu/nsimu*100),fix((acce/isimu)*100), fix((acce/(isimu+ntot))*100));
    end

    pars(1,:) = oldpar;
    scores(1) = oldss;
    for dr=1:ndelays
      pars(dr+1,:) = oldpar+randn(1,npar)*Rs{dr};     % a new proposal

      % check bounds
      if any(pars(dr+1,:)<bounds(1,:)) | any(pars(dr+1,:)>bounds(2,:))
        scores(dr+1) = Inf;

      else % inside bounds, check if accepted
        scores(dr+1)  = feval(ssfun,pars(dr+1,:),data);   % sum-of-squares

        %[alpha, Ns, Ds] = dr_acceptance(scores, pars, iRs, Ns, Ds, dr, sigma2);
        [alpha, Ns, Ds] = dr_acceptance(scores, pars, iRs, Ns, Ds, dr);

        if (log(rand) <= alpha) % we accept
          acce     = acce+1;
          oldpar   = pars(dr+1,:);
          oldss    = scores(dr+1,:);
          newi     = isimu;
          ntot     = ntot + dr - 1;
          %oldprior = newprior;
          break;
        end
      end
    end
    chain(isimu,:) = [oldpar oldss];
    % update the error variance sigma2
%    if (s20 > 0)
%      sigma2  = 1./gammar_mt(1,1,(n0+n)./2,2./(n0*s20+oldss));
%      s2chain(isimu,:) = sigma2;
%    end

    if (adaptint>0 & fix(isimu/adaptint) == isimu/adaptint)
%      beep;keyboard

      % adapt the proposal covariances
      % update covariance and mean of the chain
      [chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,1:end-1),1, ...
                                         chaincov,chainmean,wsum);
      lasti = isimu;
      [Ra,is] = chol(chaincov + eye(npar)*qcovadj);
      if (is) % singular cmat
        fprintf('Warning cmat singular, not adapting\n');
      else
        R = Ra*adascale;
        for i=1:ndelays
          Rs{i} = R./(drscale^(i-1)); % second proposal for DR try
          iRs{i} = inv(Rs{i});
        end
      end
    end

    if (do_log)
      [fid, err] = fopen(['.' filesep log_file '.dat'], 'a');
      fprintf(fid, [uuid '%ld : %e |'], isimu, oldss);
      fprintf(fid, ' %f',oldpar);
      fprintf(fid, '\n');
      fclose(fid);
    end

    if (stall_thresh > 0 & isimu-newi > stall_thresh)
      warning(['Likelihood had not evolved for ' num2str(ceil(stall_thresh)) ' iterations, aborting']);
      nsimu = isimu;
      chain = chain(1:nsimu, :);

      break;
    end
  end

  % calculate covariance and mean of the chain
  [chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,1:end-1),1, ...
                                     chaincov,chainmean,wsum);

  results.class = 'MCMC results';
  results.accepted=acce./nsimu;              % acceptance ratio
  results.mean = chainmean;
  results.cov  = chaincov;
  results.qcov = R'*R;
  results.R = R;
  results.nsimu = nsimu;
  results.drscale = drscale;
  results.adascale = adascale;
  results.adaptint = adaptint;

  return;
end

%%%%%%%%
function y=getpar(options,par,default)
%GETPAR get parameter value from a struct
% options   options struct
% par       parameter value to extract from the struct
% default   default value if par is not a member of the options struct

  if (isfield(options,par))
    y = getfield(options,par);
  elseif (nargin>2)
    y = default;
  else
    error(sprintf('Need value for option: %s',par));
  end

  return;
end

function y=gammar_mt(m,n,a,b)
%GAMMAR_MT random deviates from gamma distribution
% 
%  GAMMAR_MT(M,N,A,B) returns a M*N matrix of random deviates from the Gamma
%  distribution with shape parameter A and scale parameter B:
%
%  p(x|A,B) = B^-A/gamma(A)*x^(A-1)*exp(-x/B)
%
%  Uses method of Marsaglia and Tsang (2000)

% G. Marsaglia and W. W. Tsang:
% A Simple Method for Generating Gamma Variables,
% ACM Transactions on Mathematical Software, Vol. 26, No. 3,
% September 2000, 363-372.

% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.3 $  $Date: 2006/04/25 09:58:29 $

  if (nargin < 4)
    b=1;
  end

  y = zeros(m,n);
  for j=1:n
    for i=1:m
      y(i,j) = gammar_mt1(a,b);
    end
  end

  return;
end

function y=gammar_mt1(a,b)

  if (a<1)
    y = gammar_mt1(1+a,b)*rand(1)^(1/a);
  else
    d = a-1/3;
    c = 1/sqrt(9*d);
    while(1)
      while(1)
        x = randn(1);
        v = 1+c*x;
        if (v > 0)
          break;
        end
      end
      v = v^3;
      u = rand(1);
      if (u < 1-0.0331*x^4)
        break;
      end
      if (log(u) < 0.5*x^2+d*(1-v+log(v)))
        break;
      end
    end
    y = b*d*v;
  end

  return;
end

function [xcov,xmean,wsum]=covupd(x,w,oldcov,oldmean,oldwsum)
%COVUPD covariance update
% [xcov,xmean,wsum]=covupd(x,w,oldcov,oldmean,oldwsum)


% Marko Laine <Marko.Laine@Helsinki.FI>
% $Revision: 1.2 $  $Date: 2004/06/16 05:04:28 $

  [n,p]=size(x);
  if (n == 0) % nothing to update with
    xcov = oldcov;
    xmean = oldmean;
    wsum = oldwsum;

    return
  end

  if (nargin<2 | isempty(w))
    w = 1;
  end
  if (length(w) == 1)
    w = ones(n,1)*w;
  end

  if (nargin>2 & ~isempty(oldcov)) % update
    for i=1:n
      xi     = x(i,:);
      wsum   = w(i);
      xmeann = xi;
      xmean  = oldmean + wsum/(wsum+oldwsum)*(xmeann-oldmean);

      xcov =  oldcov + ...
              wsum./(wsum+oldwsum-1) ...
              .* (oldwsum/(wsum+oldwsum) ...
                  .* ((xi-oldmean)' *(xi-oldmean))  ...
                  - oldcov);
      wsum    = wsum+oldwsum;
      oldcov  = xcov;
      oldmean = xmean;
      oldwsum = wsum;
    end
    
  else % no update

    wsum  = sum(w);
    xmean = zeros(1,p);
    xcov  = zeros(p,p);
    for i=1:p
      xmean(i) = sum(x(:,i).*w)./wsum;
    end
    if (wsum>1)
      %%% (wsum-oldwsum/wsum)
      for i=1:p
        for j=1:i
          xcov(i,j) = (x(:,i)-xmean(i))' * ((x(:,j)-xmean(j)).*w)./(wsum-1);
          if (i ~= j)
            xcov(j,i) = xcov(i,j);
          end
        end
      end
    end

  end

  return;
end
