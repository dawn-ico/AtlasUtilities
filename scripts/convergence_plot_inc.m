function convergence_plot_inc(fname)
  clf;
  conv = csvread(fname,2,0);
  hold on;
  x = log(conv(:,1)); Linf = log(conv(:,2)); L1 = log(conv(:,3)); L2 = log(conv(:,4));  
  plot(x,Linf,'LineWidth',3);
  plot(x,L1,'LineWidth',3);
  plot(x,L2,'LineWidth',3);
  plot(x,-2*x+8,'LineWidth',3,'k-.');
  plot(x,ones(size(x))*2,'LineWidth',3,'k--');
  set(gca, 'xdir', 'reverse');
  hcb = legend(['L_{\infty}';'L_1';'L_2';'const.';'-2']);
  xlabel('log(\Delta x)[°]')
  ylabel('log(L) [-]')
  set(gca,'FontSize',50)
  set(hcb,'FontSize',50)
  ylim([-2,8]);
 endfunction
