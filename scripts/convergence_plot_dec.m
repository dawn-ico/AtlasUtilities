function convergence_plot_dec(fname)
  clf;
  conv = csvread(fname,2,0);
  hold on;
  x = log(conv(:,1)); Linf = log(conv(:,2)); L1 = log(conv(:,3)); L2 = log(conv(:,4));  
  plot(x,Linf,'LineWidth',3);
  plot(x,L1,'LineWidth',3);
  plot(x,L2,'LineWidth',3);
  plot(x,x-5.5,'LineWidth',3,'k--');
  plot(x,2*x-10,'LineWidth',3,'k-.');
  set(gca, 'xdir', 'reverse');
  hcb = legend(['L_{\infty}';'L_1';'L_2';'lin. conv';'quad. conv']);
  xlabel('log(\Delta x)[°]')
  ylabel('log(L) [-]')
  set(gca,'FontSize',20)
  set(hcb,'FontSize',20)
  ylim([-6,0]);
 endfunction
