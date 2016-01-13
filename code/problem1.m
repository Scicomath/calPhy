% 第一题 用四点高斯求积公式求二重积分

node = [0.3399810, 0.8611363];
coef = [0.6521452, 0.3478548];

u = [-node(2),-node(1),node(1),node(2)];
v = u;
A = [coef(2),coef(1),coef(1),coef(2)];

s = 0;

for i = 1:4
    for j = 1:4
	s = s + A(i)*A(j)*log(0.3*u(i) + 0.5*v(j) + 4.2);
    end
end

I = 0.075 * s;

format long
disp(I)
