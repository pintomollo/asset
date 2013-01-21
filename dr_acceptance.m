function [acc_ratio, Ns, Ds] = dr_acceptance(scores, params, iR, Ns, Ds, ntrials, n)

  %%%% Assuming gaussian distributions and integrating variance out !!!
  dist = -0.5*norm((params(ntrials+1,:) - params(ntrials,:))*iR{1})^2;
  %N = -0.5*((scores(ntrials+1)/sigma2) + dist);
  %D = -0.5*((scores(ntrials)/sigma2) + dist);
%  if (integrate)
    % Should normally be ^(-n/2) but we assume score is squared already
%    N = (log(scores(ntrials+1))*(-n/4) + dist);
%    D = (log(scores(ntrials))*(-n/4) + dist);
%  else
    N = -scores(ntrials+1) + dist;
    D = -scores(ntrials) + dist;
%  end

  Ns{ntrials, ntrials} = N;
  Ds{ntrials, ntrials} = D;

  for i=2:ntrials
    pos = ntrials - i + 1;

    dist = -0.5*norm((params(ntrials+1,:) - params(ntrials-i+1,:))*iR{i})^2;
%    prob = exp(-0.5*dist);

    Ds{ntrials, pos} = dist + Ds{ntrials-1, pos} + log(1 - min(exp(Ns{ntrials-1, pos} - Ds{ntrials-1, pos}), 1));
    Ns{ntrials, pos} = dist + Ns{ntrials, pos+1} + log(1 - min(exp(Ds{ntrials, pos+1} - Ns{ntrials, pos+1}), 1));
    %Ds{ntrials, pos} = prob*(Ds{ntrials-1, pos} - Ns{ntrials-1, pos});
    %Ns{ntrials, pos} = prob*(Ns{ntrials, pos+1} - Ds{ntrials, pos+1});
  end

  acc_ratio = Ns{ntrials, 1} - Ds{ntrials, 1};

  return;
end
