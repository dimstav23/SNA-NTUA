function result = myEgoFunction(Adj,N)

for i = 1 : N 
    node_neighbors = kneighbors(Adj,i,1) %oi geitones poy apexoyn 1 apo mena
    Data = [i, nodeneighbors]; % node indice and its neighbors indices
    
    
end
