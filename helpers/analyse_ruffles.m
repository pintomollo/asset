function analyse_ruffles(groups)

  %keyboard

  close all;

  if (nargin < 1)
    groups = {'*920*wt*.mat', '*saps1*.mat'};
  end

  dt = 10;
  event_time = [-1080 -970 -380 0 60] / dt;
  pos_pc = [1 2.1 4.3 5.4];

  for i=1:length(groups)
    ls_dir = dir(groups{i});

    res(i) = struct('pos',[],'mov',[],'depth',[],'time',[],'num',[],'area',[], 'duration', [], 'names', {{}}, 'axes', [], 'cortex', []);
    dic(i) = struct('pos',[],'mov',[],'depth',[],'time',[],'num',[],'area',[], 'duration', [], 'names', {{}}, 'axes', [], 'cortex', []);
    for j=1:length(ls_dir)
      %ls_dir(j).name;
      load(ls_dir(j).name);

      if (~isfield(mymovie.markers, 'ruffles'))
        continue;
      end

      res(i) = extract_results(mymovie.markers, res(i), j);
      res(i).names{j} = ls_dir(j).name;
      dic(i) = extract_results(mymovie.dic, dic(i), j);
    end

    res(i).pos = res(i).pos - 2*pi*(res(i).pos > 2*pi);
    dic(i).pos = dic(i).pos - 2*pi*(dic(i).pos > 2*pi);
  end

  %max_area = max(max(res(1).area),max(res(2).area));

  oks1=(res(1).depth>0);
  oks2=(res(2).depth>0);
  oks3=(dic(1).depth>0);

  nbins = 40;
  pix2um = 6.45 / 63;

  max_depth = 150 * pix2um;

  figure;boxplot(res(1).axes(:,1) * pix2um);
  figure;boxplot(res(1).axes(:,2) * pix2um);
  figure;boxplot(res(1).axes(:,1)./res(1).axes(:,2));
  figure;boxplot(res(1).cortex(:,1) * pix2um^2);
  figure;boxplot(res(1).cortex(:,2) * pix2um);
  %[bin1, bin_pos] = hist(res(1).time);
  %[bin2] = histc(res(2).time, bin_pos);
  %figure;bar(bin_pos, [bin1.' bin2]);

  % All ruffles
  %range = [max(min(res(1).time),min(res(2).time)) min(max(res(1).time), max(res(2).time))];
  %event_time = [max(min(res(1).time),min(res(2).time)) event_time min(max(res(1).time), max(res(2).time))];
  

  if (false)
    
    range = event_time([1 5]);
    range1 = oks1 & res(1).time >= range(1) & res(1).time < range(2);
    range2 = oks2 & res(2).time >= range(1) & res(2).time < range(2);
    
    ends = max([res(1).depth(range1); res(2).depth(range2)]);
    [fit1, opt1] = robustfit(res(1).depth(range1), res(1).depth(range1));
    [fit2, opt2] = robustfit(res(2).depth(range2), res(2).depth(range2));

    figure;scatter(res(1).time(range1), res(1).depth(range1), 'MarkerEdgeColor', [0.5 0.5 1]);hold on;
    scatter(res(2).time(range2), res(2).depth(range2), 'MarkerEdgeColor', [1 0.5 0.5]);

    figure;scatter(res(1).depth(range1), res(1).depth(range1), 'MarkerEdgeColor', [0.5 0.5 1]);hold on;
    plot([0 ends],[0 ends]*fit1(2) + fit1(1));
    plot([0 ends],[0 ends]*(fit1(2)+opt1.se(2)) + fit1(1) + opt1.se(1));
    plot([0 ends],[0 ends]*(fit1(2)-opt1.se(2)) + fit1(1) - opt1.se(1));
    scatter(res(2).depth(range2), res(2).depth(range2), 'MarkerEdgeColor', [1 0.5 0.5]);
    plot([0 ends],[0 ends]*fit2(2) + fit2(1), 'r');
    plot([0 ends],[0 ends]*(fit2(2)+opt2.se(2)) + fit2(1) + opt2.se(1), 'r');
    plot([0 ends],[0 ends]*(fit2(2)-opt2.se(2)) + fit2(1) - opt2.se(1), 'r');
    set(gca, 'Box', 'on')

    % Ruffling ruffles

    range = event_time(1:2);
    range1 = oks1 & res(1).time >= range(1) & res(1).time < range(2);
    range2 = oks2 & res(2).time >= range(1) & res(2).time < range(2);
    
    range = event_time([3 5]);
    range1 = range1 | oks1 & res(1).time >= range(1) & res(1).time < range(2);
    range2 = range2 | oks2 & res(2).time >= range(1) & res(2).time < range(2);

    [fit1, opt1] = robustfit(res(1).depth(range1), res(1).depth(range1));
    [fit2, opt2] = robustfit(res(2).depth(range2), res(2).depth(range2));

    figure;scatter(res(1).time(range1), res(1).depth(range1), 'MarkerEdgeColor', [0.5 0.5 1]);hold on;
    scatter(res(2).time(range2), res(2).depth(range2), 'MarkerEdgeColor', [1 0.5 0.5]);

    figure;scatter(res(1).depth(range1), res(1).depth(range1), 'MarkerEdgeColor', [0.5 0.5 1]);hold on;
    plot([0 ends],[0 ends]*fit1(2) + fit1(1));
    plot([0 ends],[0 ends]*(fit1(2)+opt1.se(2)) + fit1(1) + opt1.se(1));
    plot([0 ends],[0 ends]*(fit1(2)-opt1.se(2)) + fit1(1) - opt1.se(1));
    scatter(res(2).depth(range2), res(2).depth(range2), 'MarkerEdgeColor', [1 0.5 0.5]);
    plot([0 ends],[0 ends]*fit2(2) + fit2(1), 'r');
    plot([0 ends],[0 ends]*(fit2(2)+opt2.se(2)) + fit2(1) + opt2.se(1), 'r');
    plot([0 ends],[0 ends]*(fit2(2)-opt2.se(2)) + fit2(1) - opt2.se(1), 'r');
    set(gca, 'Box', 'on')

    % Pseudo-clevage ruffles
    %range = event_time(2:3);
    %range1 = oks1 & res(1).time >= range(1) & res(1).time < range(2);
    %range2 = oks2 & res(2).time >= range(1) & res(2).time < range(2);
    
    %[fit1, opt1] = robustfit(res(1).depth(range1), res(1).depth(range1));
    %[fit2, opt2] = robustfit(res(2).depth(range2), res(2).depth(range2));

    %figure;scatter(res(1).depth(range1), res(1).depth(range1), 'MarkerEdgeColor', [0.5 0.5 1]);hold on;
    %plot([0 ends],[0 ends]*fit1(2) + fit1(1));
    %plot([0 ends],[0 ends]*(fit1(2)+opt1.se(2)) + fit1(1) + opt1.se(1));
    %plot([0 ends],[0 ends]*(fit1(2)-opt1.se(2)) + fit1(1) - opt1.se(1));
    %scatter(res(2).depth(range2), res(2).depth(range2), 'MarkerEdgeColor', [1 0.5 0.5]);
    %plot([0 ends],[0 ends]*fit2(2) + fit2(1), 'r');
    %plot([0 ends],[0 ends]*(fit2(2)+opt2.se(2)) + fit2(1) + opt2.se(1), 'r');
    %plot([0 ends],[0 ends]*(fit2(2)-opt2.se(2)) + fit2(1) - opt2.se(1), 'r');
    %set(gca, 'Box', 'on')

    %[a, s, g] = mymean(res(1).depth(oks1), [], res(1).time(oks1));
    %figure;errorbar(g,a,s);hold on;
    %[a, s, g] = mymean(res(2).depth(oks2), [], res(2).time(oks2));
    %errorbar(g,a,s, 'r');
  end

  %keyboard

  bounds = [0:nbins]*2*pi/nbins;
  bin_pos = bounds(1:end-1) + bounds(2)/2;
  bounds(end) = bounds(end) + 1e-2;

  %figure;scatter(res(1).pos(oks1), res(1).depth(oks1));
  %hold on;scatter(res(2).pos(oks2), res(2).depth(oks2), 'r')

  range = event_time(1:2);
  range1 = oks1 & res(1).time >= range(1) & res(1).time < range(2);
  range2 = oks2 & res(2).time >= range(1) & res(2).time < range(2);

  %figure;scatter(res(1).pos(range1), res(1).length(range1), 'MarkerEdgeColor', [0.5 0.5 1]);
  %hold on;scatter(res(2).pos(range2), res(2).length(range2), 'MarkerEdgeColor', [1 0.5 0.5]);

  %keyboard
  

  figure;scatter(res(1).pos(range1), res(1).depth(range1)*pix2um, 'MarkerEdgeColor', [0.5 0.5 1]);
  hold on;scatter(res(2).pos(range2), res(2).depth(range2)*pix2um, 'MarkerEdgeColor', [1 0.5 0.5]);

  pts1 = res(1).depth(range1)*pix2um;
  [junk, assign] = histc(res(1).pos(range1), bounds);
  mean1 = zeros(nbins,1);
  m = mymean(pts1, [], assign);
  mean1(junk(1:end-1)~=0) = m;
  pts2 = res(2).depth(range2)*pix2um;
  [junk, assign] = histc(res(2).pos(range2), bounds);
  mean2 = zeros(nbins,1);
  m = mymean(pts2, [], assign);
  mean2(junk(1:end-1)~=0) = m;
  %figure;bar(bounds(1:end-1), [mean1 mean2])
  plot(bin_pos, mean1, 'LineWidth', 3);
  hold on;plot(bin_pos, mean2, 'r', 'LineWidth', 3);
  ylim([0 max_depth]);
  xlim([0 2*pi])
  set(gca, 'XTick', [0:pi/2:2*pi]);
  set(gca, 'Box', 'on')

  keyboard

  range = event_time(2:3);
  range1 = oks1 & res(1).time >= range(1) & res(1).time < range(2);
  range2 = oks2 & res(2).time >= range(1) & res(2).time < range(2);
  range3 = oks3 & dic(1).time >= range(1) & dic(1).time < range(2);

  figure;scatter(res(1).pos(range1), res(1).depth(range1), 'MarkerEdgeColor', [0.5 0.5 1]);
  hold on;scatter(dic(1).pos(range3), dic(1).depth(range3), 'MarkerEdgeColor', [0.5 1 0.5]);

  pts1 = res(1).depth(range1);
  [junk, assign] = histc(res(1).pos(range1), bounds);
  mean1 = zeros(nbins,1);
  m = mymean(pts1, [], assign);
  mean1(junk(1:end-1)~=0) = m;

  pts2 = dic(1).depth(range3);
  [junk, assign] = histc(dic(1).pos(range3), bounds);
  mean2 = zeros(nbins,1);
  m = mymean(pts2, [], assign);
  mean2(junk(1:end-1)~=0) = m;
  %figure;bar(bounds(1:end-1), [mean1 mean2])
  plot(bin_pos, mean1, 'LineWidth', 3);
  hold on;plot(bin_pos, mean2, 'g', 'LineWidth', 3);
  set(gca, 'XTick', [0:pi/2:2*pi]);
  ylim([0 max_depth]);
  xlim([0 2*pi])
  set(gca, 'Box', 'on')



  figure;scatter(res(1).pos(range1), res(1).depth(range1), 'MarkerEdgeColor', [0.5 0.5 1]);
  hold on;scatter(res(2).pos(range2), res(2).depth(range2), 'MarkerEdgeColor', [1 0.5 0.5]);
  %keyboard

  %pts1 = res(1).depth(range1);
  %[junk, assign] = histc(res(1).pos(range1), bounds);
  %mean1 = zeros(nbins,1);
  %m = mymean(pts1, [], assign);
  %mean1(junk(1:end-1)~=0) = m;
  pts2 = res(2).depth(range2);
  [junk, assign] = histc(res(2).pos(range2), bounds);
  mean2 = zeros(nbins,1);
  m = mymean(pts2, [], assign);
  mean2(junk(1:end-1)~=0) = m;
  %figure;bar(bounds(1:end-1), [mean1 mean2])
  plot(bin_pos, mean1, 'LineWidth', 3);
  hold on;plot(bin_pos, mean2, 'r', 'LineWidth', 3);
  plot(repmat(pos_pc,2,1), repmat([0 max_depth]',1,4),'k');
  set(gca, 'XTick', [0:pi/2:2*pi]);
  ylim([0 max_depth]);
  xlim([0 2*pi])
  set(gca, 'Box', 'on')

  pc_indx1 = range1 & ~((res(1).pos > pos_pc(1) & res(1).pos < pos_pc(2)) | (res(1).pos > pos_pc(3) & res(1).pos < pos_pc(4)));
  pc_indx2 = range2 & ~((res(2).pos > pos_pc(1) & res(2).pos < pos_pc(2)) | (res(2).pos > pos_pc(3) & res(2).pos < pos_pc(4)));

  invs1 = [res(1).time(pc_indx1) res(1).mov(pc_indx1)];
  [junk, junk, frame1] = unique(invs1, 'rows');
  counts1 = histc(frame1, [0.5:max(frame1)+0.5]);
  invs2 = [res(2).time(pc_indx2) res(2).mov(pc_indx2)];
  [junk, junk, frame2] = unique(invs2, 'rows');
  counts2 = histc(frame2, [0.5:max(frame2)+0.5]);

  [errs(1,1), stds(1,1)] = mymean(counts1);
  [errs(1,2), stds(1,2)] = mymean(counts2);
  [H(1), p(1)] = ttest2(counts1, counts2);
  
  %[errs(2,1), stds(2,1)] = mymean(res(1).depth(pc_indx1) ./ res(1).depth(pc_indx1));
  %[errs(2,2), stds(2,2)] = mymean(res(2).depth(pc_indx2) ./ res(2).depth(pc_indx2));
  %[H(2), p(2)] = ttest2(res(1).depth(pc_indx1) ./ res(1).depth(pc_indx1), res(2).depth(pc_indx2) ./ res(2).depth(pc_indx2));

  [errs(2,1), stds(2,1)] = mymean(res(1).depth(pc_indx1));
  [errs(2,2), stds(2,2)] = mymean(res(2).depth(pc_indx2));
  [H(2), p(2)] = ttest2(res(1).depth(pc_indx1), res(2).depth(pc_indx2));
  [errs(3,1), stds(3,1)] = mymean(res(1).area(pc_indx1));
  [errs(3,2), stds(3,2)] = mymean(res(2).area(pc_indx2));
  [H(3), p(3)] = ttest2(res(1).area(pc_indx1), res(2).area(pc_indx2));

  %length_indx1 = (res(1).length > 0 & pc_indx1);
  %length_indx2 = (res(2).length > 0 & pc_indx2);
  %[errs(5,1), stds(5,1)] = mymean(res(1).length(length_indx1));
  %[errs(5,2), stds(5,2)] = mymean(res(2).length(length_indx2));
  %[H(5), p(5)] = ttest2(res(1).length(length_indx1), res(2).length(length_indx2));

  errs
  stds
  H
  p

  %figure;scatter(res(2).mov(pc_indx2), res(2).length(pc_indx2));
  %beep;keyboard

  if (false)

  all_features = [counts1 / errs(1,1); ...
  counts2/ errs(1,1); ...
  res(1).depth(pc_indx1) ./ (res(1).depth(pc_indx1) * errs(2,1)); ...
  res(2).depth(pc_indx2) ./ (res(2).depth(pc_indx2) * errs(2,1)); ...
  res(1).depth(pc_indx1)/ errs(3,1); ...
  res(2).depth(pc_indx2)/ errs(3,1); ...
  res(1).depth(pc_indx1)/ errs(4,1); ...
  res(2).depth(pc_indx2)/ errs(4,1)];
  %res(1).length(length_indx1)/ errs(5,1); ...
  %res(2).length(length_indx2)/ errs(5,1)];

  groups = [1*ones(size(counts1)); ...
  2*ones(size(counts2)); ...
  3*ones(size(res(1).depth(pc_indx1) ./ res(1).depth(pc_indx1))); ...
  4*ones(size(res(2).depth(pc_indx2) ./ res(2).depth(pc_indx2))); ...
  5*ones(size(res(1).depth(pc_indx1))); ...
  6*ones(size(res(2).depth(pc_indx2))); ...
  7*ones(size(res(1).depth(pc_indx1))); ...
  8*ones(size(res(2).depth(pc_indx2)))];
  %9*ones(size(res(1).length(length_indx1))); ...
  %10*ones(size(res(2).length(length_indx2)))];

  figure;
  boxplot(all_features, groups);
  end

  stds = stds ./ errs(:,[1 1])
  errs = errs ./ errs(:,[1 1])
  figure;
  barweb(errs, stds, [], {'Number', 'Ratio', 'Depth', 'Area', 'Length'}, '', 'Property', 'Ratio with WT',[],'y',[],1);


  range = event_time(3:4);
  range1 = oks1 & res(1).time >= range(1) & res(1).time < range(2);
  range2 = oks2 & res(2).time >= range(1) & res(2).time < range(2);

  figure;scatter(res(1).pos(range1), res(1).depth(range1), 'MarkerEdgeColor', [0.5 0.5 1]);
  hold on;scatter(res(2).pos(range2), res(2).depth(range2), 'MarkerEdgeColor', [1 0.5 0.5]);

  pts1 = res(1).depth(range1);
  [junk, assign] = histc(res(1).pos(range1), bounds);
  mean1 = zeros(nbins,1);
  m = mymean(pts1, [], assign);
  mean1(junk(1:end-1)~=0) = m;
  pts2 = res(2).depth(range2);
  [junk, assign] = histc(res(2).pos(range2), bounds);
  mean2 = zeros(nbins,1);
  m = mymean(pts2, [], assign);
  mean2(junk(1:end-1)~=0) = m;
  %figure;bar(bounds(1:end-1), [mean1 mean2])
  plot(bin_pos, mean1, 'LineWidth', 3);
  hold on;plot(bin_pos, mean2, 'r', 'LineWidth', 3);
  set(gca, 'XTick', [0:pi/2:2*pi]);
  ylim([0 max_depth]);
  xlim([0 2*pi])
  set(gca, 'Box', 'on')


  range = event_time(4:5);
  range1 = oks1 & res(1).time >= range(1) & res(1).time < range(2);
  range2 = oks2 & res(2).time >= range(1) & res(2).time < range(2);

  figure;scatter(res(1).pos(range1), res(1).depth(range1), 'MarkerEdgeColor', [0.5 0.5 1]);
  hold on;scatter(res(2).pos(range2), res(2).depth(range2), 'MarkerEdgeColor', [1 0.5 0.5]);

  pts1 = res(1).depth(range1);
  [junk, assign] = histc(res(1).pos(range1), bounds);
  mean1 = zeros(nbins,1);
  m = mymean(pts1, [], assign);
  mean1(junk(1:end-1)~=0) = m;
  pts2 = res(2).depth(range2);
  [junk, assign] = histc(res(2).pos(range2), bounds);
  mean2 = zeros(nbins,1);
  m = mymean(pts2, [], assign);
  mean2(junk(1:end-1)~=0) = m;
  %figure;bar(bounds(1:end-1), [mean1 mean2])
  plot(bin_pos, mean1, 'LineWidth', 3);
  hold on;plot(bin_pos, mean2, 'r', 'LineWidth', 3);
  set(gca, 'XTick', [0:pi/2:2*pi]);
  ylim([0 max_depth]);
  xlim([0 2*pi])
  set(gca, 'Box', 'on')

  %keyboard
  event_time / 6

  return;
end

function res = extract_results(mystruct, res, indx)

  ruffles = mystruct.ruffles;
  if (isfield(mystruct,'cytokinesis'))
    zero = mystruct.cytokinesis;
    times = [1:length(ruffles)];
    cytokinesis = (times - zero);
  else
    'No cyto ?'
    return;
    %data.cytokinesis = frame_timing(mystruct);
  end

  %keyboard
  timing = [];
  links = [];

  for k=1:length(ruffles)
    pos = carth2elliptic(ruffles(k).carth, mystruct.centers(:,k), mystruct.axes_length(:,k),mystruct.orientations(:,k));

    if (isempty(pos))
      continue;
    end

    pos = pos(:,1);
    %if (data(i,j).invert)
    %  pos = pos + pi;
    %end
    
    res.pos = [res.pos; pos];
    res.mov = [res.mov; ones(size(pos))*indx];
    res.depth = [res.depth; ruffles(k).properties(:,1) + ruffles(k).properties(:,3)];
    res.area = [res.area; ruffles(k).properties(:,2)];
    %res.area = [res.area; ruffles(k).properties(:,1) + ruffles(k).properties(:,3)];
    %res.depth = [res.depth; ruffles(k).properties(:,2)];
    res.time = [res.time; ones(size(pos))*cytokinesis(k)];
    res.num = [res.num; length(pos) cytokinesis(k)];

    res.axes = [res.axes; mystruct.axes_length(:,k)'];
    res.cortex = [res.cortex; [polyarea(mystruct.cortex(k).carth(:,1), mystruct.cortex(k).carth(:,2)), sum(sqrt(sum(diff(mystruct.cortex(k).carth).^2,2))) + sum(ruffles(k).properties(:,3))]];

    new_links = ruffles(end-k+1).cluster(:,1);
    tmp = ones(size(new_links));
    tmp(links(links~=0)) = timing(links~=0) + 1; 
    timing = tmp;
    links = new_links;
    if (any(links == 0))
      res.duration = [res.duration; [timing(links==0) ones(size(links(links==0)))*cytokinesis(end-k+1)]];
    end
  end

  return;
end
