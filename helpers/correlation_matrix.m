function [rel_C, C, rel_H, H] = correlation_matrix(pts, precision)

  if (nargin < 2)
    precision = 2e-6;
  end

  legend = {'D_A', 'k_{A+}', 'k_{A-}', 'k_{AP}', '\alpha','\rho_A','\psi', 'L', ...
            'D_P', 'k_{P+}', 'k_{P-}', 'k_{PA}', '\beta', '\rho_P', '\nu', '\gamma'};

  if (isstruct(pts))

    nrates = size(pts.rate, 2);
    noffset = size(pts.offset, 2);
    nenergy = size(pts.energy, 2);
    nvisc = size(pts.viscosity, 2);
    nflow = size(pts.flow, 2);
    nscale = size(pts.flow_scaling, 2);

    pts = [pts.score pts.rate pts.offset pts.energy pts.viscosity pts.flow pts.flow_scaling];

    legend = [legend(1:nrates), ... 
              cellstr([repmat('\delta_', noffset, 1), num2str([1:noffset].')]).', ...
              cellstr([repmat('E_', nenergy, 1), num2str([1:nenergy].')]).', ...
              cellstr([repmat('\eta_', nvisc, 1), num2str([1:nvisc].')]).', ...
              legend(end-(nflow+nscale-1):end)];

    legend = legend(~cellfun('isempty', legend));

  else
    nshifts = size(pts,2) - length(legend) - 1;
    shifts = cellstr([repmat('\delta_', nshifts, 1), num2str([1:nshifts].')]);
    legend = [legend(1:end-1) shifts.' legend(end)];
  end

  is_fixed = any(isnan(pts), 1) | all(bsxfun(@eq, pts, pts(1,:)), 1);
  pts = pts(:, ~is_fixed);
  legend = legend(~is_fixed(2:end));

  scores = pts(:,1);
  pts = pts(:,2:end);
  std_values = median(pts);

  norm_values = std_values;
  norm_values(norm_values == 0) = 1;

  nparams = size(pts, 2);
  H = NaN(nparams);

  rel_H = H;
  C = H;
  rel_C = H;

  for i=1:nparams
    valsi = unique(pts(:,i));
    %dptsi = abs(valsi - std_values(i));
    nis = floor((numel(valsi)-1)/2);
    ci = nis + 1;

    for ni = 1:nis
      for j=i:nparams
        if (isnan(H(i,j)))
          pattern = std_values;
          if (pattern(i)==0 && nis>1)
          %if (false)
            %ind_1 = ci;
            %ind0 = ci+ni;
            %ind1 = ci+ni+1;
            ind_1 = ci;
            ind0 = ci+ceil((nis-ni+1)/2);
            ind1 = ci+(nis-ni+1);
          else
            %ind_1 = ci-ni;
            %ind0 = ci;
            %ind1 = ci+ni;
            ind_1 = ni;
            ind0 = ci;
            ind1 = ci+(nis-ni+1);
          end

          if (i==j)
            pattern(i) = valsi(ind_1);
            f_1 = find_value(pts, pattern, scores);

            if (isempty(f_1))
              continue;
            end

            pattern(i) = valsi(ind0);
            f0 = find_value(pts, pattern, scores);

            if (isempty(f0))
              continue;
            end

            pattern(i) = valsi(ind1);
            f1 = find_value(pts, pattern, scores);

            if (isempty(f1))
              continue;
            end

            H(i,j) = (f_1 - 2*f0 + f1) / (abs(valsi(ind0)-valsi(ind_1)) * ...
                                          abs(valsi(ind1)-valsi(ind0)));

            rel_H(i,j) = H(i,j)*norm_values(i)*norm_values(j);

            if (abs(rel_H(i,j)) < precision)
              H(i,j) = NaN;
              rel_H(i,j) = NaN;
            else
              H(j,i) = H(i,j);
              rel_H(j,i) = rel_H(i,j);
            end
          else
            valsj = unique(pts(:,j));
            %dptsj = abs(valsj - std_values(j));

            njs = floor((numel(valsj)-1)/2);
            cj = njs + 1;

            for nj = 1:njs
              if (isnan(H(i,j)))
                if (pattern(j)==0 && (njs>1))
                %if (false)
                  %jnd_1 = cj;
                  %jnd0 = cj+nj;
                  %jnd1 = cj+nj+1;
                  jnd_1 = cj;
                  jnd0 = cj+ceil((njs-nj+1)/2);
                  jnd1 = cj+(njs-nj+1);
                else
                  %jnd_1 = cj-nj;
                  %jnd0 = cj;
                  %jnd1 = cj+nj;
                  jnd_1 = nj;
                  jnd0 = cj;
                  jnd1 = cj+(njs-nj+1);
                end

                pattern(i) = valsi(ind_1);
                pattern(j) = valsj(jnd_1);
                f_1_1 = find_value(pts, pattern, scores);

                if (isempty(f_1_1))
                  continue;
                end

                pattern(i) = valsi(ind_1);
                pattern(j) = valsj(jnd1);
                f_11 = find_value(pts, pattern, scores);

                if (isempty(f_11))
                  continue;
                end

                pattern(i) = valsi(ind1);
                pattern(j) = valsj(jnd_1);
                f1_1 = find_value(pts, pattern, scores);

                if (isempty(f1_1))
                  continue;
                end

                pattern(i) = valsi(ind1);
                pattern(j) = valsj(jnd1);
                f11 = find_value(pts, pattern, scores);

                if (isempty(f11))
                  continue;
                end

                H(i,j) = (f_1_1 + f11 - f_11 - f1_1) / ...
                    (abs(valsi(ind0)-valsi(ind_1))*abs(valsj(jnd0)-valsj(jnd_1)) + ...
                     abs(valsi(ind0)-valsi(ind_1))*abs(valsj(jnd1)-valsj(jnd0)) + ...
                     abs(valsi(ind1)-valsi(ind0))*abs(valsj(jnd0)-valsj(jnd_1)) + ...
                     abs(valsi(ind1)-valsi(ind0))*abs(valsj(jnd1)-valsj(jnd0)));

                rel_H(i,j) = H(i,j)*norm_values(i)*norm_values(j);

                if (abs(rel_H(i,j)) < precision)
                  H(i,j) = NaN;
                  rel_H(i,j) = NaN;
                else
                  H(j,i) = H(i,j);
                  rel_H(j,i) = rel_H(i,j);
                end
              end
            end
          end
        end
      end
    end
  end

  %keyboard

  H(isnan(H)) = 0;
  rel_H(isnan(rel_H)) = 0;

  goods = ~all(H==0);

  %Cov = inv(H(goods, goods));
  %[eig_vect, eig_vals] = eig(Cov);

  %rel_Cov = inv(rel_H(goods, goods));
  %[rel_eig_vect, rel_eig_vals] = eig(rel_Cov);
  %keyboard;

  Cov = pinv(H(goods, goods));
  d = diag(Cov);
  s = sign(d);
  D = pinv(diag(sqrt(abs(d))));
  %correct = diag(s);
  %correct(correct==0) = 1;
  Corr = D*Cov*D;
  %Corr = Corr .* correct;
  Corr(s<0,:) = NaN;
  Corr(:,s<0) = NaN;
  C(goods, goods) = Corr;

  Cov = pinv(rel_H(goods, goods));
  d = diag(Cov);
  s = sign(d);
  D = pinv(diag(sqrt(abs(d))));
  %correct = diag(s);
  %correct(correct==0) = 1;
  Corr = D*Cov*D;
  %Corr = Corr .* correct;
  Corr(s<0,:) = NaN;
  Corr(:,s<0) = NaN;
  rel_C(goods, goods) = Corr;

  return;
end

function val = find_value(matrix, pattern, score)

  diffs = sum(bsxfun(@minus, matrix, pattern).^2, 2);
  [diff_val, indx]=min(diffs);
  %good = (sum(bsxfun(@minus, matrix, pattern).^2, 2) < 2e-6);

  if (diff_val < 2e-6)
    val = score(indx);
    %val = val(1);
  else
    val = [];
  end

  return;
end
