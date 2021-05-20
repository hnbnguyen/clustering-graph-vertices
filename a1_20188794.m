function set12 = a1_20188794(elist)
% Function for CISC271, Winter 2021, Assignment #1
%
% IN:
%     elist - Mx2 array of edges, each row is a pair of vertices
% OUT:
%     set12 - Nx1 vertex clsutering, -1 for SET1 and +1 for SET2

    % Problem size: number of vertices in the graph
    n = max(elist(:));

   % Creating adjacency matrix
A = zeros(n);
k = 4
[el_row, el_col] = size(elist);
for el_ind_row = 1: el_row
    a = elist(el_ind_row,1);
    b = elist(el_ind_row,2);
    if a ~= k & b ~= k
        
    A(a,b) = 1;
    end
end
A = zeros(n);

[e_row_num, e_col_num] = size(elist);
k = 5
for e_ind_row = 1: el_row

    if a ~= k && b ~= k

        A(E(e_ind_row,1),E(e_ind_row,2)) = 1;

    else

        pass;

    end

end
%%Create Lapacian matrix
%Creating the sum vector of all edges per row
B = sum(A,2);
%Adding the sum matrix into the adjancency matrix
for ind_diag = 1:n
    A(ind_diag, ind_diag) = B(ind_diag,1);
end
%Changing existing ones into negative ones
for row_ind = 1:n 
    for col_ind = 1:n
        if A(row_ind, col_ind) == 1
            if row_ind ~= col_ind
                A(row_ind, col_ind) = -1;
            end
        end
    end
end
disp(A)
%Finding the second vector of matrix A
[Emat, lraw] = eig(A,'vector');
v2 = Emat(:,2);
for count= 1 : n
    if v2(count,1) < 0
        v2(count,1) = -1;
    else 
        v2(count,1) = 1;
    end
end
set12 = v2;

%Sorting which vertices belong to Set 1 or Set 2
S1 = zeros(1,1);
S1_count = 1;
S2 = zeros(1,1);
S2_count = 1;
for ind_fiedler = 1:n
    if set12(ind_fiedler,1) == -1
        S1(1,S1_count) = ind_fiedler;
        S1_count = S1_count + 1;
    else
        S2(1, S2_count) = ind_fiedler;
        S2_count = S2_count + 1;
    end
end

disp('Set 1 vertices are:');
disp(S1);
disp('Set 2 vertices are:');
disp(S2);

    % Plot the graph, Cartesian and clustered
    plot271a1(A, set12);
end

function plot271a1(Amat, cvec)
% PLOTCLUSTER(AMAT,CVEC) plots the adjacency matrix AMAT twice;
% first, as a Cartesian grid, and seconnd, by using binary clusters
% in CVEC to plot the graph of AMAT based on two circles
%
% INPUTS: 
%         Amat - NxN adjacency matrix, symmetric with binary entries
%         cvec - Nx1 vector of class labels, having 2 distinct values
% OUTPUTS:
%         none
% SIDE EFFECTS:
%         Plots into the current figure

    % %
    % % Part 1 of 2: plot the graph as a rectangle
    % %

    % Problem size
    [m n] = size(Amat);

    % Factor the size into primes and use the largest as the X size
    nfact = factor(n);
    nx = nfact(end);
    ny = round(n/nx);

    % Create a grid and pull apart into coordinates; offset Y by +2
    [gx, gy] = meshgrid((1:nx) - round(nx/2), (1:ny) + 2);

    % Offset the odd rows to diagram the connections a little better
    for ix=1:2:ny
        gx(ix, :) = gx(ix, :) + 0.25*ix;
    end

    % The plot function needs simple vectors to create the graph
    x = gx(:);
    y = flipud(gy(:));

    % Plot the graph of A using the Cartesian grid
    plot(graph(tril(Amat, -1), 'lower'), 'XData', x, 'YData', y);
    axis('equal');

    % %
    % % Part 2 of 2: plot the graph as pair of circles
    % %
    % Set up the X and Y coordinates of each graph vertex
    xy = zeros(2, numel(cvec));

    % Number of cluster to process
    kset = unique(cvec);
    nk = numel(kset);

    % Base circle is radius 2, points are centers of clusters
    bxy = 2*circlen(nk);

    % Process each cluster
    for ix = 1:nk
        jx = cvec==kset(ix);
        ni = sum(jx);
        xy(:, jx) = bxy(:, ix) + circlen(ni);
    end

    hold on;
    plot(graph(Amat), 'XData', xy(1,:), 'YData', xy(2,:));
    hold off;
    title(sprintf('Clusters of (%d,%d) nodes', ...
        sum(cvec==kset(1)), sum(cvec==kset(2))));
end

function xy = circlen(n)
% XY=CIRCLEN(N) finds N 2D points on a unit circle
%
% INPUTS:
%         N  - positive integer, number of points
% OUTPUTS:
%         XY - 2xN array, each column is a 2D point

    xy = [cos(2*pi*(0:(n-1))/n) ; sin(2*pi*(0:(n-1))/n)];
end
