function analyse_centrosomes(fname)

  close all;

  if (nargin == 0)
    fname = 'GZ376_wt_MERGE_';
  end

  if (ischar(fname))
    load(fname);
  else
    mymovie = fname;
  end

  rot_thresh = [0.02 0.02];
  var_thresh = 0.2;

  % Already rescaled by ASSET !!
  %pix2um = 6.45 / 63;
  pix2um = 1;
  frame_rate = 10;
  centrosomes = [];
  avg = [];

  cytok = mymovie.dic.cytokinesis;

  full_range = 1:length(mymovie.data.centrosomes);
  for i=full_range
    centrosomes = cat(4,centrosomes, mymovie.data.centrosomes(i).all);
    avg = cat(4,avg, mymovie.data.centrosomes(i).warped);
  end
  nmov = size(centrosomes, 3);

  full_data = squeeze(any(any(squeeze(sum(isnan(centrosomes), 3)) <= 1, 1),2));
  first = find(full_data, 1, 'first');
  last = find(full_data, 1, 'last');

  range = [first:last];
  plot_range = (range - cytok) / 6;

  centr_avg = diff(avg(:,:,:,:), [], 1);
  centr_avg = squeeze(atan2(-centr_avg(1,2,:,:), centr_avg(1,1,:,:)));
  centr_avg = centr_avg + 2*pi*(centr_avg < 0);
  centr_avg = fix_angles(centr_avg);

  rot_acc = abs(diff(centr_avg, 2));

  %meeting = find(abs(centr_avg - pi/2) < angl_thresh, 1, 'last');
  %displac = find((abs(centr_avg - pi) < angl_thresh) & (full_range' > first), 1, 'first');

  %keyboard

  % Display angular rotation
  centr_angl = diff(centrosomes, [], 1);
  centr_angl = squeeze(atan2(-centr_angl(1,2,:,:), centr_angl(1,1,:,:)));
  centr_angl = centr_angl + 2*pi*(centr_angl < 0);
  centr_angl = fix_angles(centr_angl);
  [junk, stds] = mymean(centr_angl, [], 1);
  var_ok = (stds(2:end-1) < var_thresh);
  var_ok = var_ok(2:end-1).';

  rot_avg = mean([rot_acc(1:end-2) rot_acc(2:end-1) rot_acc(3:end)],2);

  %figure;
  %plot(plot_range,rot_acc(first:last))
  %hold on;plot(plot_range,stds,'g');
  %keyboard

  %meeting = find(rot_acc < rot_thresh, 1, 'first') + 1;
  %displac = find(rot_acc < rot_thresh, 1, 'last') + 1;
  meeting = find(rot_avg < rot_thresh(1) & var_ok, 1, 'first') + 2;
  displac = find(rot_avg < rot_thresh(2) & var_ok, 1, 'last') + 2;

  el1 = squeeze(avg(1,1,1,:));
  el2 = squeeze(avg(2,1,1,:));

  %[junk, nebd] = min(abs(diff(el1)) + abs(diff(el2)) + 100*(full_range(1:end-1)' < first));

  centering_range = [meeting:displac];

  %centr_avg = centr_avg(range);
  %keyboard

  mymap = ((range - 1) / range(end-1)).';
  mymap = repmat(mymap, [1, 2, 3]);
  myred = repmat(reshape([1 0 0],[1 1 3]),[length(range) 2 1]);
  myblue = repmat(reshape([0 0 1],[1 1 3]),[length(range) 2 1]);
  myred2 = repmat(reshape([1 0.5 0.5],[1 1 3]),[length(range) 2 1]);
  myblue2 = repmat(reshape([0.5 0.5 1],[1 1 3]),[length(range) 2 1]);

  display_reference;hold on;
  for i=1:nmov
    color_line(squeeze(centrosomes(1,1,i,range)),-squeeze(centrosomes(1,2,i,range)),mymap.*myblue2, 'LineWidth', 1);
    color_line(squeeze(centrosomes(2,1,i,range)),-squeeze(centrosomes(2,2,i,range)),mymap.*myred2, 'LineWidth', 1);
  end
  color_line(squeeze(avg(1,1,1,range)),-squeeze(avg(1,2,1,range)),mymap.*myblue, 'LineWidth', 3);
  color_line(squeeze(avg(2,1,1,range)),-squeeze(avg(2,2,1,range)),mymap.*myred, 'LineWidth', 3);

  % Display angular rotation
  figure;hold on;
  plot([0 3*pi/2; 0 3*pi/2]', [plot_range(meeting-first) plot_range(meeting-first); plot_range(displac-first) plot_range(displac-first)]','k')
  %keyboard
  for i=1:nmov
    plot(centr_angl(i,range), plot_range, 'Color', [0.5 0.5 0.5]);
  end

  [junk, stds] = mymean(centr_angl(:, range), [], 1);

  plot(centr_avg(range)+stds', plot_range, 'r', 'LineWidth', 3);
  plot(centr_avg(range)-stds', plot_range, 'r', 'LineWidth', 3);
  plot(centr_avg(range), plot_range, 'k', 'LineWidth', 3);
  set(gca, 'YDir', 'reverse');
  set(gca, 'XTick', [0:pi/2:3*pi/2], 'XGrid', 'on');
  set(gca, 'Box', 'on');
  ylim([plot_range(1) plot_range(end)]);
  xlim([0 3*pi/2]);
  
  % Display centrosomes distance
  figure;hold on;
  plot([0 25; 0 25]', [plot_range(meeting-first) plot_range(meeting-first); plot_range(displac-first) plot_range(displac-first)]','k')
  tmp = zeros(nmov, length(range));
  for i=1:nmov
    dist = diff(centrosomes(:,:,i,range));
    dist = squeeze(hypot(dist(:,2,:,:), dist(:,1,:,:)));
    plot(dist, plot_range, 'Color', [0.5 0.5 0.5]);

    tmp(i,:) = dist';
  end

  [junk, stds] = mymean(tmp, 1);
    
  dist = diff(avg(:,1:2,:,:));
  dist = squeeze(hypot(dist(:,2,:,:), dist(:,1,:,:)));

  plot(dist(range)+stds', plot_range, 'g', 'LineWidth', 2);
  plot(dist(range)-stds', plot_range, 'g', 'LineWidth', 2);
  plot(dist(range), plot_range, 'k', 'LineWidth', 2);
  ylim([plot_range(1) plot_range(end)]);
  xlim([0 max(dist(:))]);
  set(gca, 'YDir', 'reverse');
  set(gca, 'Box', 'on')

  % Display in AP distance
  figure;hold on;
  plot([0 0], [0 plot_range(end)], 'k');
  plot([-25 25; -25 25]', [plot_range(meeting-first) plot_range(meeting-first); plot_range(displac-first) plot_range(displac-first)]','k')
  for i=1:nmov
    %cel1 = 100 * (1 - ((squeeze(centrosomes(1,1,i,range)) + 25) / 50));
    %cel2 = 100 * (1 - ((squeeze(centrosomes(2,1,i,range)) + 25) / 50));
    cel1 = squeeze(centrosomes(1,1,i,range));
    cel2 = squeeze(centrosomes(2,1,i,range));

    plot(cel1, plot_range, 'Color', [0.5 0.5 1]);
    plot(cel2, plot_range, 'Color', [1 0.5 0.5]);
  end

  [junk, stds] = mymean(centrosomes(:,1,:,range),3);

  plot(el1(range)+stds(1,:)', plot_range, 'g', 'LineWidth', 2);
  plot(el1(range)-stds(1,:)', plot_range, 'g', 'LineWidth', 2);
  plot(el2(range)+stds(2,:)', plot_range, 'g', 'LineWidth', 2);
  plot(el2(range)-stds(2,:)', plot_range, 'g', 'LineWidth', 2);
  plot(el1(range), plot_range, 'b', 'LineWidth', 2);
  plot(el2(range), plot_range, 'r', 'LineWidth', 2);
  xlim([-25 25])
  ylim([plot_range(1) plot_range(end)]);
  %set(gca, 'XDir', 'reverse');

  set(gca, 'YDir', 'reverse');
  set(gca, 'Box', 'on')

  % Rescale to um/min movements
  diff1 = squeeze(diff(avg(1,:,1,:),[],4)) / frame_rate;
  diff2 = squeeze(diff(avg(2,:,1,:),[],4)) / frame_rate;

  urange = 0.5;

  figure;
  for i=1:nmov
    cdiff1 = squeeze(diff(centrosomes(1,:,i,:),[],4)) / frame_rate;
    cdiff2 = squeeze(diff(centrosomes(2,:,i,:),[],4)) / frame_rate;

    subplot(1,3,1);hold on;
    if (i==1)
      plot([-urange urange], [0 0], 'k');
      plot([0 0], [-urange urange], 'k');
    end
    scatter(cdiff1(1,first:meeting-1),-cdiff1(2,first:meeting-1),'MarkerEdgeColor',[0.5 0.5 1]);
    scatter(cdiff2(1,first:meeting-1),-cdiff2(2,first:meeting-1),'MarkerEdgeColor',[1 0.5 0.5]);
    subplot(1,3,2);hold on;
    if (i==1)
      plot([-urange urange], [0 0], 'k');
      plot([0 0], [-urange urange], 'k');
    end
    scatter(cdiff1(1,meeting:displac-1),-cdiff1(2,meeting:displac-1),'MarkerEdgeColor',[0.5 0.5 1]);
    scatter(cdiff2(1,meeting:displac-1),-cdiff2(2,meeting:displac-1),'MarkerEdgeColor',[1 0.5 0.5]);
    subplot(1,3,3);hold on;
    if (i==1)
      plot([-urange urange], [0 0], 'k');
      plot([0 0], [-urange urange], 'k');
    end
    scatter(cdiff1(1,displac:last-1),-cdiff1(2,displac:last-1),'MarkerEdgeColor',[0.5 0.5 1]);
    scatter(cdiff2(1,displac:last-1),-cdiff2(2,displac:last-1),'MarkerEdgeColor',[1 0.5 0.5]);
  end

  subplot(1,3,1);hold on;
  scatter(diff1(1,first:meeting-1),-diff1(2,first:meeting-1),'b', 'LineWidth', 2);
  scatter(diff2(1,first:meeting-1),-diff2(2,first:meeting-1),'r', 'LineWidth', 2);
  axis equal
  xlim([-urange urange])
  ylim([-urange urange])
  set(gca, 'Box', 'on')

  subplot(1,3,2);hold on;
  scatter(diff1(1,meeting:displac-1),-diff1(2,meeting:displac-1),'b', 'LineWidth', 2);
  scatter(diff2(1,meeting:displac-1),-diff2(2,meeting:displac-1),'r', 'LineWidth', 2);
  axis equal
  xlim([-urange urange])
  ylim([-urange urange])
  set(gca, 'Box', 'on')

  subplot(1,3,3); hold on;
  scatter(diff1(1,displac:last-1),-diff1(2,displac:last-1),'b', 'LineWidth', 2);
  scatter(diff2(1,displac:last-1),-diff2(2,displac:last-1),'r', 'LineWidth', 2);
  axis equal
  xlim([-urange urange])
  ylim([-urange urange])
  set(gca, 'Box', 'on')

  %centr_mov = diff(centrosomes, [],4) * pix2um;
  %centr_mov = squeeze(hypot(centr_mov(:,1,:,:), centr_mov(:,2,:,:)));

  avg_mov = diff(avg, [],4) * 60 / frame_rate;
  avg_mov = squeeze(hypot(avg_mov(:,1,:,:), avg_mov(:,2,:,:)));
  %avg_mov = mean(avg_mov)

  avg_rot = diff(centr_avg) * 60 / frame_rate;
  centr_rot = diff(centr_angl, [], 2);

  tmp = avg_mov(:,first:meeting-1);
  [me(1,1:2), sd(1,1:2)] = mymean(tmp, 2);
  tmp = avg_mov(:,meeting:displac-1);
  [me(2,1:2), sd(2,1:2)] = mymean(tmp, 2);
  tmp = avg_mov(:,displac:last-1);
  [me(3,1:2), sd(3,1:2)] = mymean(tmp, 2);

  tmp = avg_rot(first:meeting-1);
  [me(1,3), sd(1,3)] = mymean(tmp(:));
  tmp = avg_rot(meeting:displac-1);
  [me(2,3), sd(2,3)] = mymean(tmp(:));
  tmp = avg_rot(displac:last-1);
  [me(3,3), sd(3,3)] = mymean(tmp(:));

  dist = (diff(dist) * 60 / frame_rate)';
  tmp = dist(first:meeting-1);
  [me(1,4), sd(1,4)] = mymean(tmp, 2);
  tmp = dist(meeting:displac-1);
  [me(2,4), sd(2,4)] = mymean(tmp, 2);
  tmp = dist(displac:last-1);
  [me(3,4), sd(3,4)] = mymean(tmp, 2);


  me
  sd
  %figure;
  %barweb(me.', sd.', [], {'Movements', 'Rotation'});

  [splindle, stds] = mymean(abs(diff(centrosomes(:,1,:,last))),3)

  %size(centrosomes)
  %size(avg)
  avg_diff = centrosomes - repmat(avg(:,1:2,:,:), [1 1 nmov 1]);
  avg_err = squeeze(hypot(avg_diff(:,2,:,first:last),avg_diff(:,2,:,first:last)));
  [err, stds] = mymean(avg_err(:))

  ([first meeting displac last] - cytok) / 6

  return;
end

function new_angl = fix_angles(angl)

  %keyboard

  [nmov, npts] = size(angl);
  new_angl = angl;
  tests = repmat([0 2*pi -2*pi], nmov, 1);

  for i=npts-1:-1:1
    %if (i==39)
    %  keyboard
    %end

    new_pts = angl(:,[i i i]) + tests;
    dist = abs(new_pts - new_angl(:,[i+1 i+1 i+1]));
    [tmp, indxs] = min(dist, [], 2);
    %keyboard
    new_angl(:, i) = new_pts(sub2ind([nmov 3],[1:nmov]',indxs));
  end

  return;
end
