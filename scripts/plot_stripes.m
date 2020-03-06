function plot_stripes
  hold on;
  for i=0:31
    fname = sprintf('stripe_%04d.txt',i);
    stripe = load(fname);
    scatter3(stripe(:,1),stripe(:,2),stripe(:,3),50,'filled');
  endfor
end