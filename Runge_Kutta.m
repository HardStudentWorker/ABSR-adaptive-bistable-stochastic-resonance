function x = Runge_Kutta(a,b,h,s)
%使用四阶Runge-Kutta方程求解随机共振非线性系统的输出
%输入为双稳态系统的参数a&b、迭代步长h
%输出为x
x = zeros(length(s),1);
x(1) = -sqrt(a/b);

for i=1:length(s)-1
    k1 = h * (a*x(i) - b*x(i).^3 + s(i));
    k2 = h * (a*(x(i) + k1/2) - b*(x(i) + k1/2).^3 + (s(i)+s(i+1))/2);%这里做了一定修正，原为s(i)
    k3 = h * (a*(x(i) + k2/2) - b*(x(i) + k2/2).^3 + (s(i)+s(i+1))/2);%这里做了一定修正，原为s(i+1)
    k4 = h * (a*(x(i) + k3) - b*(x(i) + k3).^3 + s(i+1));
    x(i+1) = x(i) + (k1 + 2*k2 + 2*k3 +k4)/6;
end

end