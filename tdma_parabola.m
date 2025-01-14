dx = 0.001;
dt = dx^2;
alp = dt / (dx^2);
k = floor(1 / dx);
u = zeros(1, k + 1);
time = 0;

tic;  
while time < 0.01 %0.02,0.05,0.1,0.2
    
    tdma = zeros(k + 1, k + 1);
    tdma(1, 1) = 1;
    for i = 2:k
        tdma(i, i - 1) = -alp;
        tdma(i, i) = 1 + 2 * alp;
        tdma(i, i + 1) = -alp;
    end
    tdma(k + 1, k) = -2 * alp;
    tdma(k + 1, k + 1) = 1 + 2 * alp;

   
    r = u;

   
    for i = 2:k + 1
        f = tdma(i, i - 1) / tdma(i - 1, i - 1);
        r(i) = r(i) - f * r(i - 1);
        
        if i == k + 1
            tdma(k + 1, k + 1) = tdma(k + 1, k + 1) - f * tdma(k, k + 1);
            continue;
        end
        for j = i - 1:i + 1
            tdma(i, j) = tdma(i, j) - f * tdma(i - 1, j);
        end
    end

   
    u(k + 1) = r(k + 1) / tdma(k + 1, k + 1);
    u(1) = 1;  
    for i = k:-1:2
        u(i) = (r(i) - tdma(i, i + 1) * u(i + 1)) / tdma(i, i);
    end

    
    time = time + dt;
    disp(time);
end
toc; 


plot(u);
xlabel('Position Index');
ylabel('U Value');
title('U vs X');


