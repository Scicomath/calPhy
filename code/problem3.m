% 第三题 用有限差分法求解单位园域中的泊松方程
% 设置初始值
omega = 1.25;
epsilon = 10^(-7);
r0 = 1.0;
h = 0.1;
deltaR = h;
global N = r0/h+1;
global M = 16;
deltaPhi = 2*pi/(M-1);

R = 0:deltaR:r0;
Phi = 0:deltaPhi:2*pi;
U = zeros(N,M);
% 定义下标i的转换函数, 依据周期性条件
function newi = subi(i)
  global N;
  if i==0
    newi = N-1;
  elseif i == N+1
    newi = 2;
  elseif i == N
    newi = 1;
  else
    newi = i;
  end
endfunction
% 定义下标j的转换函数, 依据周期性条件
function newj = subj(j)
  global M;
  if j==0
    newj = M-1;
  elseif j == M+1
    newj = 2;
  elseif j == M
    newj = 1;
  else
    newj = j;
  end
endfunction

function y = fun(r, phi)
  y = -50 * r^2 * sin(2*phi);
endfunction
% 进行迭代
do
  maxDiffU = 0;
  for i = (N-1):-1:1 
    if i==1 % 判断是否为圆心
      oldU = U(1,1);
      meanU = mean(U(2,1:(M-1)));
      U0 = omega*(meanU - h^2/4*fun(R(i),Phi(j))) + (1-omega)*oldU;
      U(1,:) = U0;
      diffU = abs(U(1,1)-oldU);
      if diffU > maxDiffU
	maxDiffU = diffU;
      end
    else 
      for j = 1:M-1
	oldU = U(i,j);
	alpha_0 = (deltaPhi*(i-1))^(-2);
	alpha_1 = 1 - (2*i-2)^(-1);
	alpha_2 = 1 + (2*i-2)^(-1);
	fenmu = 2*(1+alpha_0);

	temp1 = alpha_0 * ( U(subi(i),subj(j-1)) + \
			    U(subi(i),subj(j+1)) );
	temp2 = alpha_1 * U(subi(i-1),subj(j));
	temp3 = alpha_2 * U(subi(i+1),subj(j));
	U(i,j) = omega*(temp1 + temp2 + temp3 - \
			h^2*fun(R(i),Phi(j)) )/fenmu + (1-omega)*U(i,j);
	diffU = abs(U(i,j) - oldU);
	if diffU > maxDiffU
	  maxDiffU = diffU;
	end
      end
    U(i,M) = U(i,1);
    end
  end
until maxDiffU < epsilon
% 计算解析解
realU = zeros(size(U));
for i = 1:N
  for j = 1:M
    realU(i,j) = 25/6*(R(i))^2*(1-(R(i))^2)*sin(2*Phi(j));
  end
end
% 输出结果
disp(max(max(abs(U-realU))))
save -ascii U.dat U
save -ascii realU.dat realU
