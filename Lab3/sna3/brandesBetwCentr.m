% Brandes Algorithm for Betweenness Centrality %
% Author: Konstantinos Sotiropoulos %
% Import: NXN adjacency matrix %
% Export: 1XN betweenness centrality for every graph node



function BC=brandesBetwCentr(adjMatrix)

nodesNumber = length( adjMatrix );
%% Pre-Processing of Adjacency Matrix
pSize = max(sum(adjMatrix));
testMatrix = zeros(nodesNumber, pSize);
nodesNumber = length( adjMatrix );
indexMatrix = zeros(nodesNumber,1);
 for i=1:nodesNumber
     tempArray = find(adjMatrix(i,:));
     for j=1:length( tempArray )
         testMatrix(i,j) = tempArray(j);
     end
     indexMatrix(i) = length( tempArray );
     
 end
clear adjMatrix tempArray
P = zeros( nodesNumber,pSize );
tic
%% initialize betw, begin actual code
BC = zeros( nodesNumber, 1 );

%distance = ones( nodesNumber, 1).*(-1);
% for s c- V do
for source=1:nodesNumber
    % S<- empty stack;
    Stack = zeros( 1, nodesNumber );
    stackIndex = nodesNumber;
    
    % P[w]<- empty list, w c- V;
    indexP = ones( nodesNumber,1 );
    
    sigma = zeros( nodesNumber, 1);
    sigma( source ) =1;
    
    distance = ones( nodesNumber, 1).*(-1);
    distance( source ) = 0;
    
    % Q is a queue
    qStart = 1;
    qEnd = 1;
    Q =  [source] ;
    while ( qStart<=qEnd )
       
        v = Q( qStart );
        qStart = qStart+1;
        
        Stack(stackIndex) = v;
        stackIndex = stackIndex-1;
        
        neighboursV = testMatrix(v,1:indexMatrix(v));
        for n=1:indexMatrix(v) 
            w = neighboursV(n);
            % w found for the first time?
            if distance(w)<0
               qEnd = qEnd + 1;
               Q( qEnd ) = w;
               distance(w) = distance(v)+1;
            end
            % shortest path to w via v?
            if ( distance(w) ==  distance(v) + 1 )
                sigma(w) = sigma(w)+sigma(v);
                P(w,indexP(w)) = v;
                indexP(w) = indexP(w)+1;
            end
        end
    end
    delta = zeros( nodesNumber, 1 );
    % S returns vertices in order of non-increasing distance from s
    for node=1:nodesNumber
        % pop w<-S
        w = Stack( node );
        indexP(w) = indexP(w)-1;
        for j=1:indexP(w)
            v = P(w,j);
            delta(v) = delta(v)+( sigma(v)/sigma(w) )*(1 + delta(w)) ;
        end
        
        if ( w~=source )
             BC(w) = BC(w)+delta(w);
        end
    end
    
end
% if graph is undirected divide by 2
BC = BC / 2;
toc
clear Stack Q distance w v source distance
