function CS= fitness_function(Sub,r)
	%Calculates the Community Score of subgraph S of our graph 
	us = sum(sum(Sub));
	
	M = sum(mean(Sub,2).^r)/length(Sub);
	
	CS = M*us;
end
