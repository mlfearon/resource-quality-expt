function sumall = calc_r_Fx(r, Fx, lx, x)

x=x(1:end-1);
lx=lx(1:end-1);
erx=exp(-r*x);

sumall = sum(Fx.*lx.*erx)-1;
sumall = abs(sumall);

