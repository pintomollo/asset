function weight = edges_intens(img, params)

  weight = 1 - imadm_mex(img);
  weight(isnan(weight)) = Inf;

  return;
end
