function sumall = calc_r_mx(r, mx, lx, x)

erx=exp(-r*x);

sumall = sum(lx.*mx.*erx)-1;
sumall = abs(sumall);

