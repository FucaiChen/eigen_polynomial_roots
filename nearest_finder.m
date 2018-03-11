% newton-down hill test
clear ; close all; clc
% polycoeff = [ 0     1       0.1          -0.2    10;...
%                         0      0.7   -0.2       0.03     0.001;...
%                         0,     0       0.0        0.00    0.0];
% pos = [1, 0, 0]';
polycoeff=[   0, 1, -2, 3;
						0, -0.5, 1, 2;
						0, 0, 0, 0]
pos = [2, 2, 2]';
order = size(polycoeff, 2) - 1;


t = -0.4: 0.01: 1;
for i = 1: (order + 1)
        T(i, :) = t.^(i - 1);
end
p = polycoeff * T;

d = p - [pos(1) * ones(1, length(t)); pos(2) * ones(1, length(t)); pos(3) * ones(1, length(t)); ];
D = d(1,:).^2 + d(2,:).^2 + d(3,:).^2 ;
D_diff = (D(2:end) - D(1:end - 1) )/ (t(2) - t(1));

figure(1)
plot(p(1,:), p(2,:));
for i =int32( linspace(1,  length(t), 10));
    text(p(1,i),p(2,i) ,num2str(t(i)));
end
hold on; grid on; axis equal;
plot(pos(1), pos(2), '+r')
figure
plot(t, D);
hold on; grid on;
plot(t(2:end), D_diff)
legend('D', 'D_diff')
clear t p T d D D_dot
%%

range = [0, 1];
init_t =0.1;
t_k = init_t;
C = polycoeff;
p = pos;
Q = C; 
Q(:, 1) = Q(:, 1) - p;
J = Q' * Q;
T              = zeros(order + 1, 1);
dTdt       = zeros(order + 1, 1);
d2Tdt2   = zeros(order + 1, 1);
dDdt = 0;
d2Ddt2 = 0;



n = order + 1;
eig_value = 0: (n-1);
m_temp = diag(eig_value);
M = J * m_temp;
M_temp = M;
M_temp(:,1) = [];
p = diag_sum(M_temp);
c =  fliplr(p);
while (1)
    lambda = 1;
    t = t_k;
    for i = 1: (order + 1)
        T(i)            = t^(i - 1);
        dTdt(i)      = (i - 1) * t^max(i - 2, 0);
        d2Tdt2(i) = max(i - 2, 0) * (i - 1) * t^max(i - 3, 0);
    end
    dDdt = dTdt' * J * T + T' * J *dTdt;
    d2Ddt2 = d2Tdt2' * J * T + dTdt' * J * dTdt + dTdt' * J * dTdt + T' * J * d2Tdt2;
    f_tk = dDdt;
    f_tk_dot = d2Ddt2;
    
    %% test if f_tkp1 > f_tk
    while(1)
        t_kp1 = t_k - lambda * f_tk /f_tk_dot;
        t = t_kp1;
        for i = 1: (order + 1)
            T(i)            = t^(i - 1);
            dTdt(i)      = (i - 1) * t^max(i - 2, 0);
            d2Tdt2(i) = max(i - 2, 0) * (i - 1) * t^max(i - 3, 0);
        end
        dDdt = dTdt' * J * T + T' * J *dTdt;
        d2Ddt2 = d2Tdt2' * J * T + dTdt' * J * dTdt + dTdt' * J * dTdt + T' * J * d2Tdt2;
        f_tkp1 = dDdt;
        if abs(f_tkp1) < abs(f_tk)
            break;
        else
            lambda = 0.5 * lambda;        
        end   
    end
     if abs(t_kp1 - t_k) < 0.000001
        break;
     end
     t_k = t_kp1
    end
        t_k


        
        
        
    
