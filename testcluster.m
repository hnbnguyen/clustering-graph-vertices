%function testcluster(n)

% Create an adjacency matrix with 2 equally sized clusters,
% having a "few" connections between them
%Araw = adjgen(n, .4*n, .8*n);

% Randomly permute the nodes of the graph
B = adjperm(Araw);

% Find the edge list of A
[jx,ix] = find(tril(B));
elist = [ix jx];

% TEST CODE: try knocking out some list entries
if 0
    checkval = 2;
    elist(ismember(elist(:,1), checkval+[0]),:) = [];
    elist(ismember(elist(:,2), checkval+[0]),:) = [];
end

% Find the adjacency matrix of the edge list
A = zeros(size(B));
A(elist(:,1) + (elist(:,2) - 1)*size(A,1)) = 1;
A = A + A';

% Compute the degree matrix and Laplacian matrix of B
D = diag(sum(A, 2));
L = D - A;

% Cluster the graph into 2 sets: Fiedler vector
[evecs, evals] = eig(L, 'vector');
kF = (kmeans(evecs(:,2),2) > 1.5)*2 - 1;

% Plot the graph, Cartesian and clustered
plot271a1(A, kF);
