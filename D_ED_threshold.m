function threshold = D_ED_threshold(sigma_n,N,PF,p_BSR)
%此函数用于仿真确定检测阈值
%输入：
% sigma_n   估计得到的噪声参数
% N         检测点数
% PF        设定的虚警率
% p_BSR     最优BSR系统的参数[a,b,h]
%输出：
% threshold D-ED检测算法的阈值

%-----------流程------------
%由噪声参数生成M个长度为N的噪声片段-->送入给定参数的BSR系统中计算
% 输出的检验统计量-->对检验统计量从大到小排序,第M*PF个即为阈值

%噪声序列的个数
M = 10000;

%生成噪声
% noise = sigma_n * randn(M,N);%每一行为一个噪声序列

%通过BSR系统输出,并分别计算D-ED检验统计量
a = p_BSR(1);
b = p_BSR(2);
h = p_BSR(3);
T_y = zeros(M,1);
parfor i = 1:M
    y = Runge_Kutta(a,b,h,sigma_n*randn(N,1));
    y_mean = mean(y);
    T_y(i) = sum((y-y_mean).^2);
end

%对检验统计量从大到小排序
T_y = sort(T_y,'descend');

%输出满足虚警率的阈值
threshold = T_y(M*PF);

end