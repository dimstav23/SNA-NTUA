function [ communities ] = my_ga(Adj,pc,pm,elitism,generations,chrom,fit_max_unchanged,order_for_fitness_function)
	%my genetic algorithm for communities detection
	%%INPUT
	%Adj : adjacency matrix
	%pc : crossover probability
	%pm : mutation probability
	%elitism : number of conquering first chromosomes for the next generation
	%generations : number of generations
	%chrom : number of chromosomes
	%fit_max_unchanged : how many times fitness function maximum remains the same
	
	%%OUTPUT
	%communities : communities as a cell array
	
	%%Initialization
	% n = numer of nodes
	n=size(Adj,1);
	B = cell(chrom,generations+1);
	%get the neighbours of each node and initialize our first generation
	%taking randomly one of them 
	neighbours = cell(n);
	for i=1:n
		neighbours{i} = find(Adj(i,:)==1);
	end
	t = 1; %t represents our generation
	for i=1:chrom
		for j=1:n
			B{i,t}(j) = neighbours{j}(randi(size(neighbours{j},2)));
		end
	end
	
	%%Main Body of my Genetic Algorithm
	
	%create the new adjacency matrix according to the chromosomes to find the connected components	
	Adj_2 = zeros(n,n);
	%number of times max of fitness functions remains unchanged
	max_fitness_cnt = 0;
	max_fitness = -inf;
	while (t <= generations+1 && max_fitness_cnt < fit_max_unchanged) 
		
		for i = 1:chrom
			for j = 1:n
				Adj_2(j,B{i,t}(j)) = 1;
				Adj_2(B{i,t}(j),j) = 1;
			end
			%we now find the connected components
			conn_comp = find_conn_comp(Adj_2);
			%calculation of fitness function and scores
			r=1;
			for j = 1:size(conn_comp,2)
				%take the appropriate subgraph
				Sub{j} = subgraph(Adj_2,conn_comp{j});
				%calculation of the CS of this component
				CStemp(j) = fitness_function(Sub{j},order_for_fitness_function);
			end
			CS(i) = sum(CStemp);
			%clearvars CStemp  %%%%%%%%%%%%%%%%%%%%%%%%%%
			Adj_2 = zeros(n,n); %%%%%%%%%%%%%%%%%%%%%%%%%%%%
		end
		
		%Selection with elitism
		[sorted_CS,indexes] = sort(CS,'descend');
		%keep the  #elitism first chromosomes 
		for i=1:elitism
			B{i,t+1} = B{indexes(i),t};
		end
		%roulette method
		total_CS= sum(CS);
		for i=elitism+1:chrom %continue with the remaining chromosomes
			x=rand(1);
			k=1;
			while (k < chrom & x < sum(CS(1:i))/total_CS )
				k=k+1;		
			end
			B{i,t+1} = B{k,t}; %should i take same chromosomes or not??????? %%%%%%%%%%%%%%%%%%%%%%%
		end
		
		%one-point Crossover
		for i=1:2:chrom-1
			if (rand(1) <= pc)
				pos = randi(n-1);
				for k=pos+1:n
					aux = B{i,t+1}(k);
					B{i,t+1}(k) = B{i+1,t+1}(k);
					B{i+1,t+1}(k) = aux;				
				end
			end
		end
		  
		%Mutation
		for i=1:chrom
			for k=1:n
				if (rand(1) <= pm)
					B{i,t+1}(k) = neighbours{k}(randi(size(neighbours{k},2)));
				end
			end
		end
		
		if (sorted_CS(1) == max_fitness)
			max_fitness_cnt = max_fitness_cnt + 1;
		else
			max_fitness = sorted_CS(1);
			max_fitness_cnt = 1;
		end
		
		%break;
		%clearvars CS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		t = t+1
	end
	%get my results
	t = t - 1;
	Adj_2 = zeros(n,n);
	for j = 1:n
		Adj_2(j,B{indexes(1),t}(j)) = 1;
		Adj_2(B{indexes(1),t}(j),j) = 1;
	end
	communities = find_conn_comp(Adj_2);
end