function [delay, ref] = synchronize_domains(reference, domain)

  %keyboard

  if (size(reference,2) == 1)
    position_reference = [1:length(reference)].';

    [ref, indexes] = unique(reference, 'first');
    position_reference = position_reference(indexes);

    if (ref(1)~=0 && position_reference(1) > 1)
      ref = [0; ref];
      position_reference = [position_reference(1)-1; position_reference];
    elseif (ref(1)==0)
      position_reference(1) = find(ref==0, 1, 'last');
    end

    goods = isfinite(ref);
    ref = ref(goods);
    position_reference = position_reference(goods);
  else
    position_reference = reference(:,1);
    ref = reference(:,2);
  end

  domain(isnan(domain)) = 0;
  position_domain = [1:length(domain)].';
  goods = false(size(domain));

  prev_indx = find(isfinite(domain), 1, 'first');
  init = prev_indx + 1;
  for i=init:length(domain)
    if (domain(i) > domain(prev_indx))
      goods(i) = true;
      prev_indx = i;
    end
  end

  domain = domain(goods);
  position_domain = position_domain(goods);

  if (domain(1)~=0 && position_domain(1) > 1)
    domain = [0; domain];
    position_domain = [position_domain(1)-1; position_domain];
  elseif (domain(1)==0)
    position_domain(1) = find(domain==0, 1, 'last');
  end

  pos_dom = interp1(domain, position_domain, ref);

  %figure;hist(pos_dom - position_reference, 20);figure;
  diff_pos = pos_dom - position_reference;
  [mval, sval] = mymean(diff_pos);
  delay = median(diff_pos(diff_pos <= mval + 3*sval & diff_pos >= mval - 3*sval));

  if (isempty(delay))
    delay = 0;
  end

%  tmp1 = mean(pos_dom - position_reference);
%  tmp2 = [pos_dom ones(size(pos_dom))] \ position_reference;
%  tmp3 = robustfit(pos_dom, position_reference);
%  tmp4 = pos_dom(1) - position_reference(1);
%  tmp5 = median(pos_dom - position_reference);

%  delay = tmp3(1);

%  delay = -[delay, tmp1, tmp2(2), tmp3(1) tmp4 tmp5];

%  keyboard

  ref = [position_reference ref];

  return;
end
