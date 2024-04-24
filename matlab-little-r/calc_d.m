function NLL = calc_d(d, death)

% AUTHOR: 	Spencer Hall, sprhall@indiana.edu
% DATE: 	20 Feb 2012
% PROJECT: 	life table work -- multiple projects
% PURPOSE: 	used to estimate d, instantaneous death rate, using life table data

% INPUTS:	d is the guess used by the optimizer for d
%		death is a matrix -- structure:		
%			column 1 has time information
%			column 2 has 0 for "still alive" and 1 for "died during life table"

% OUTPUTS:	negative log likelihood, NLL

% BASED ON:	Mccallum's parameter estimation book.



NLL = 0;					% initialize the NLL				



for j = 1:1:length(death)			% no, loop through the data

   if death(j,2) == 0				% if the animal did live,...
      p = exp(-d*death(j,1));			% ... use the survival function -- see Mccallum

   else						% if the animal died during the experiment
      p = d*exp(-d*death(j,1));			% ... use the death function -- see McCallum

   end

   NLL = NLL + log(p);				% add up the log-likelihood (LL)

end 						% j loop

if d < 0					% if d is a negative number, make it a really big NLL
   NLL = 10^6;
else
   NLL = -NLL;					% else, return the negative LL
end