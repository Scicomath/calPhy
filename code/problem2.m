% 第二题 应用龙格-库塔法求初值问题
% 设置初始值
t0 = 0;
h = 0.1;

y1_0 = 1;
y2_0 = -1;

function K = fun(t, y)
  f1 = y(2);
  %f2 = exp(t) - 2*t*y(2) - t^2*y(1);
  f2 = exp(2*t)*sin(t) - 2*y(1) + 2*y(2);
  K = [f1,f2];
endfunction

t = 0:h:1;
n = size(t,2);
Y = zeros(n,2);
Y(1,1) = y1_0;
Y(1,2) = y2_0;

for i = 1:n-1
  K1 = fun(t(i), Y(i,:));
  K2 = fun(t(i) + h/2, Y(i,:) + (h/2)*K1);
  K3 = fun(t(i) + h/2, Y(i,:) + (h/2)*K2);
  K4 = fun(t(i) + h, Y(i,:) + h*K3);
  Y(i+1,:) = Y(i,:) + (K1 + 2*K2 + 2*K3 + K4) * (h/6);
end
% 画图
plot(t,Y(:,1)','b',t,Y(:,2)','g')
xlabel('t')
ylabel('y')
legend('y-t',"y'-t")
title("y-t and y'-t")


