% 第四题 边界元法求解泊松方程的混合边值问题
% 初始化变量
n = 3;
H = zeros(n);
G = zeros(n);
% 计算相关参数
function [d, rd, theta, r1, r2, s] = parafun(i,j)
endPointx = [0,1; 1,0; 0,0];
endPointy = [0,0; 0,1; 1,0];
midPointx = [0.5, 0.5, 0];
midPointy = [0, 0.5, 0.5];
xi = midPointx(i);
yi = midPointy(i);
x1 = endPointx(j,1);
x2 = endPointx(j,2);
y1 = endPointy(j,1);
y2 = endPointy(j,2);
r1Vec = [x1-xi, y1-yi];
r2Vec = [x2-xi, y2-yi];
r21Vec = [x2-x1, y2-y1];
s = sqrt(sum(r21Vec.^2));
r1 = sqrt(sum(r1Vec.^2));
r2 = sqrt(sum(r2Vec.^2));
d = -dot(r1Vec, r21Vec./s);
nVec = [(y2-y1)/s, (x1-x2)/s];
rd = dot(r1Vec, nVec);
theta = acos(dot(r1Vec,r2Vec)/(r1*r2));
endfunction

sVec = [1, sqrt(2), 1];

for i = 1:n
  for j = 1:n
    if i == j % 计算对角元素
      H(i,i) = -pi;
      G(i,i) = sVec(i)*(log(sVec(i)/2) - 1);
    else % 计算非对角元素
      [d, rd, theta, r1, r2, s] = parafun(i,j);
      H(i,j) = theta;
      G(i,j) = (s-d)*log(r2) + d*log(r1) - s + abs(rd)*theta;
    end
  end
end

% 计算A和R
A = [G(:,1:2),-H(:,3)];
R = [H(:,1:2),-G(:,3)]*[0;1;0];
uq = inv(A)*R;
disp(uq)
