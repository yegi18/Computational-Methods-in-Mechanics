u = zeros(102, 102);
de = 0.01;


u(:, 1) = 1;
u(1, :) = 1;

for i = 1:101
    u(i, 102) = -1 * de + u(i, 100);
end

for i = 1:101
    u(102, i) = (0.5 + 0.25 * (i-1)*de) * (2 * de) + u(100, i);
end

do = zeros(102, 102);
maxdo = 1;
iter = 0;


while maxdo >= 1e-6 && iter <= 100000
    v=u;
   
    for i = 2:101
        for j = 2:101
            do(i, j) = (v(i, j-1) + v(i, j+1) + v(i-1, j) + v(i+1, j) - 4 * v(i, j));
            u(i, j) = u(i, j) + do(i, j)*(1/4);
        end
    end

  
    for i=1:101
        u(102,i)= (0.5 + 0.25 * (i-1)*de) * (2 * de) + u(100, i);
    end


    for i = 1:101
        u(i, 102) = -1 * de + u(i, 100);
    end
 
    maxdo = max(max(abs(do)));
    iter = iter + 1;
    fprintf('%f, %d\n', maxdo, iter);
end


figure;
contourf(u(1:101, 1:101), 20, 'LineColor', 'none'); 
colormap(parula); 
colorbar; 
xlabel('X-axis');
ylabel('Y-axis');
title('Contour Plot of 2D Matrix');
%ITERATIONS=54582