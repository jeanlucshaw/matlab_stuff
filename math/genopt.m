function MOM = genopt(Q,X,dom,NP,NI,pc,pm)
% Synthax:         MOM = genopt(@Q,@X,dom,NP,NI,pc,pm)
%
% Optimisation routine following a genetic algorithm. Returns the best parameter
% set which solves the problem expressed by the anonymous function 'Q', when 
% either the criteria expressed by the anonymous function 'X' has been met or 
% NI iterations have been completed. 'X' should be formulated in terms of the
% values returned by 'Q'. The 'dom' variable is a matrix of size [p 2], where
% p is the number of parameters passed to the function 'Q'. It's rows are the 
% domains from which the initial population will be generated for each parameter.
% Variable 'NP' is the size of the population and must be an even number. Variables
% 'pc' and 'pm' are respectively the probability of mixing genetic material and
% probability of mutation. Recommended values are 0.8 and 0.2 .

%INIT POPULATION
NA	 = numel(dom(:,1)) ;
PP	 = nan(NP,NA) ;
SC	 = nan(NP,1) ;
RK	 = nan(NP,1) ;
for ii	 = 1:NP
	for jj = 1:NA
		PP(ii,jj) = dom(jj,1) + (dom(jj,2)-dom(jj,1))*rand(1) ;
	end
end

for kk = 1:NI 
	
	%EVALUATION
	for ii = 1:NP ; SC(ii) = Q(PP(ii,:)) ; end 
	
	%RANKING
	[HS,RK]	= sort(SC) ;
	PP	= PP(RK,:) ;
	
	%TEST
	if X(HS(1)) ; break; end
	
	%SELECTION
	DICE    = randi([2 4]) ;
	MOM 	= PP(1,   :) ;
	DAD	= PP(DICE,:) ;

	%REPRODUCTION
	PP(1,:) = MOM ;
	PP(2,:) = DAD ;
	for ii = 3:2:NP
		for jj = 1:NA
			% CROISEMENT
			if rand(1) < pc
				DICE		 = rand(1) ;
				PP(ii,jj)	 = MOM(jj)*DICE + DAD(jj)*(1-DICE) ;
				PP(ii+1,jj)	 = DAD(jj)*DICE + MOM(jj)*(1-DICE) ;
			else
				PP(ii,jj)	 = MOM(jj) ;
				PP(ii+1,jj)	 = DAD(jj) ;
			end
		
			% MUTATION
			if rand(1) < pm 
				PP(ii,jj) = dom(jj,1) + (dom(jj,2)-dom(jj,1))*rand(1) ;
			end
		end
	end
end

end
