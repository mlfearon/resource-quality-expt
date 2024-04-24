function clone_calc_rbd

% AUTHOR: 	Spencer Hall, sprhall@indiana.edu
% DATE: 	16 Feb 2012
% PROJECT: 	life table analysis work, several manuscripts
% PURPOSE: 	calculate little r two ways, little d, then b by their addition (for both r metrics)
%		follows McCallum's Estimating Population Parameters book

% 		THIS IS THE MASTER FUNCTION
%		This function inputs the data, does the calculation for each clone, then saves the result	


% INPUTS:	A text file -- change the name in the code 
%		Top row -- first cell is zero, then the x values
% 		Rest of matrix -- the life table data (zeros and babies, 99 for day of death and afterwards)
% 			First column is clone ID
%			rest of the columns = data


% OUTPUT:	A text file -- change the name in the code





% (1). Input data from a text file, preliminary work with it, and other preliminary stuff


load -ascii ResourceQuality_Metsch_LifeTable_clone_x_infstatus.txt;			% CHANGE FILE NAME HERE!!
data = ResourceQuality_Metsch_LifeTable_clone_x_infstatus;

n=data(end,1);						% last clone number; assumes 1 through n to analyse
x=data(1,2:end);					% the days of the life table (x in life table speak)

BS = 1000;						% number of bootstraps for SE estimates





% (2). Now, scroll through each clone


for j = 1:1:n						% scroll from clone 1 to n

   j							% output clone ID for diagnosis purposes


   % (A) get data for each clone

   ind = find(data(:,1) == j);				% ... get the data in two steps
   clone = data(ind, 2:end);				% ... now "clone" has the fecundity data for the focal clone
   [r c] = size(clone);					% ... and size of clone data, for bootstrapping



   % (B) point estimate of parameters

   res = calc_rbd(clone, x);				% use "calc_rbd.m" to calculate statistics for that clone
 


   % (C) now, bootstrap some SEs
 
   if BS > 1  						% calculate bootstraps as long as BS isn't <<1>>

      for k = 1:1:BS
         order=round((r-1)*rand(r,1)+1);		% (1) create order in which to take residuals
         rnd_clone=clone(order(:),:);		
         rnd_res(k,:) = calc_rbd(rnd_clone, x);
      end						% end r loop through bootstraps

      SE = std(rnd_res);

   else							% just do this if bootstrapped not wanted
      SE = [];
   end


   results(j,:) = [j, res, SE];				% package the results


end




% (3) Store output as a text file

save ResourceQuality_Metsch_lifetable_clone_x_infstatus_results.txt results -ascii		% store output   