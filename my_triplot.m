function my_triplot(fname);
  fnameP = sprintf('%sP.txt', fname);
  fnameT = sprintf('%sT.txt', fname);
  P = load(fnameP);
  T = load(fnameT);
  triplot(T,P(:,1),P(:,2));
endfunction
