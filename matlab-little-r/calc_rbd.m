function res = calc_rbd(data, x)

% AUTHOR: 	Spencer Hall, sprhall@indiana.edu
% DATE: 	16 Feb 2012
% PROJECT: 	life table analysis work, several manuscripts
% PURPOSE: 	calculate little r two ways, little d, then b by their addition (for both r metrics)
%		follows McCallum's Estimating Population Parameters book	


% INPUTS:	For a given clone

%		"data" has beaker (single animal) in a row; columns are fecundity
%		fecundities are zeros or number of babies that day
%		"99" is entered for the day of death and all following days -- the code needs that
% 		for a particular application with Stu's data, I have entered "88" for days after last censored observation

% 		"x" has the vector of ages of the animal

% OUTPUTS	res contains two birth rate estimates, the death rate estimate, and the two little rs
%		assumes that death rate, d, is zero for clones in which all animals survived experiment
%		final entry is arithmetic fecunidty, averaged for clone

%		so order is (packaged in a vector)
%			instantaneous b (mx method)	b_mx 
%			instantaneous b (Fx method)	b_Fx 			
%			little r, (mx method)		r_mx 
%			little r, (Fx method)		r_Fx 
%			instanteaneous d		d
%			mean arithmetic fecunidty	Mean_AF 
%			mean age at first reproduction	Mean_AFR






% (1). Preliminary information on the data matrix

[r,c]=size(data);



% (2). PREPARE FOR DEATH RATE CALCULATIONS
% calcualte time of death for each individual, if they died


nodie = 0;					% counter for number of animals dead

for n=1:1:r					% scroll through rows (beakers)

   death(n,1) = 99;				% dummy for death counts below
   death(n,2) = 0;

   for k = 1:1:c				% scroll through columns (x values = days)

      if data(n,k) == 99			% if death that day
         if data(n,k-1) ~= 99			% ... but not the day before
	    death(n,1) = x(k);			% ... record the day of death
            death(n,2) = 1;
         end

      elseif data(n,k) == 88			% if no die but stopped recording after certain number of clutches (coded by '88')
         if data(n,k-1) ~= 88			% ... but only do it on the first time '88' comes up
            death(n,1) = x(k);			% ... record the day of last clutch monitored (last time censused still alive)
            death(n,2) = 0;			% ... and indicate the animal is still alive
            nodie = nodie + 1;			% add another one to the 'still alive' category
         end

      end 					% if loop

   end 						% k loop

   if death(n,1) == 99				% now, if no one died and the data aren't coded with the '88' or they were monitored to the last day of expt
      death(n,1) = x(end);			% ... but the last day of the experiment in there
      nodie = nodie + 1;   			% ... and count how many animals didn't die
   end	

end 						% n loop


if nodie == r					% if no animals died		
   need_d = 0;					% ... tell code below to not calculate d, don't need to
else 
   need_d = 1;					% ... but if at least animal died, calculate d
end







% (3). SURVIVAL AND FECUNDITY SCHEDULES
% calculate lx (proportion survived) and mx (average number of offspring per day that survived)

for k = 1:1:c					% scroll through number of columns

   nalive=0; fec = 0;				% initializations

   for n = 1:1:r
      if (data(n,k) ~= 99)
         nalive = nalive + 1;			% number of hosts alive

         if data(n,k) >= 100
            indfec = data(n,k) - 100;		% I add 100 to data when animal births on day of death
         else
            if data(n,k) ~= 88
               indfec = data(n,k);
            else
               indfec = 0;
            end
         end

         fec = fec + indfec;			% add up fecundity for the day

      end 					% if statement  
   end    					% n loop

   lx(k) = nalive/r;				% fraction of initial hosts still alive

   if nalive == 0				% fecundity of alive animals
      mx(k) = 0;				% ... zero if none are alive to prevent NaN
   else
      mx(k) = fec/nalive;   
   end						% if loops

end						% end of k loop








% (3B). CALCULATE ALTERNATIVE FECUNDITY METRIC, FX
% now calculate px (probability of surviving between classes) and Fx (the possibly improved fecundity metric)

for k = 1:1:c-1
   if (lx(k+1) == 0) & (lx(k) == 0)			% if there are no animals anymore
      px(k) = 0;					% ... call px zero to avoid NaN
   else
      px(k) = lx(k+1)/lx(k);				% probability of survival from x to x+1
   end
   Fx(k) = (mx(k) + px(k)*mx(k+1))/2;			% Caswell's Fx, from McCallum pg 141
end



% (3C). CALCULATE AVERAGE FECUNDITY SCALED ARITHMETICALLY AND AGE AT FIRST REPRODUCTION



for n = 1:1:r						% scroll through beakers
   babies = 0;						% initialize baby counter
   AFRtemp = 100;					% put in a really big number for the temporary holder of the age at first repro

   for k = 1:1:c					% now scroll through days

      if (data(n,k) ~= 99) & (data(n,k) ~= 88)		% ... and add up number of babies
         babies = babies + data(n,k);
         
         if data(n,k) > 0				% if the animals reproduced
            clutchday = x(k);				% ... get the date of reproduction
            if clutchday < AFRtemp			% ... and if this is the first date of reproduction
               AFRtemp = clutchday;			% ... record this as the age at first reproduction
            end
         end
      end
   end							% end k loop though columns (days)

   ArithFec(n) = babies/death(n);				% calculate arithmetic fecunidty for each clone
   AFR(n) = AFRtemp;
end							% end n loop (through beakers)







% (4). DO THE ACUTAL r, b, and d calculations
% now, estimate little r, using mx and then Fx

r0 = 0.2; d0 = 0.02; 					% preliminary parmeter guesses



% little r calcualtions

r_mx = fminsearch('calc_r_mx', r0, [], mx, lx, x);	% the mx method
r_Fx = fminsearch('calc_r_Fx', r0, [], Fx, lx, x);	% the Fx method




% little d calculations (death rate)

if need_d == 0						% if no one died in the life table for the clone
   d = 0;						% ... set death rate equal to zero
else
   d = fminsearch('calc_d', d0, [], death);		% ... otherwise, estimate exponential death rate
end





% little b calculations (instantaneous birth rate) and mean arithmetic fecundity and mean age at first reproduction (AFR)

b_mx = r_mx + d;					% using the mx method
b_Fx = r_Fx + d;					% using the Fx method


Mean_AF = mean(ArithFec);				% mean arithmetic fecundity

ind = find(AFR ~= 100);
Mean_AFR = mean(AFR(ind));






% package the results into "res" 

res = [b_mx b_Fx r_mx r_Fx d Mean_AF Mean_AFR];