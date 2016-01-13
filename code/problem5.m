% 第五题 采用蒙特卡罗法解方程

a = 1;
% 定义函数
function y = fun1(x)
  y = exp(-x^3) - tan(x) + 800;
endfunction

function y = fun2(x)
  y = x + 5*exp(-x) - 5;
endfunction

function root = findRoot(x0, b, N, epsilon, fun)
  m = 0;
  x = x0;

  f = feval(fun,x);
  if abs(f) < epsilon
    return
  end
  
  do
    xtemp = x + b*(2*rand()-1);
    ftemp = feval(fun,xtemp);
    m = m + 1;
    if m == N
      b = b/2;
      m = 0;
    end
    if abs(ftemp) < abs(f)
      f = ftemp;
      x = xtemp;
    end
  until abs(f)<epsilon
  root = x;
endfunction

f = @fun1;
root1 = findRoot(1, 0.2, 20, 10^(-5), f);
disp(root1)
disp(fun1(root1))

f = @fun2;
root2 = findRoot(5, 0.3, 20, 10^(-5), f);
disp(root2)
disp(fun2(root2))
