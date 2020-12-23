function x_ml = J8_coding(x, Constellation_points, Nt, snr, type)
% 输入：x：发射的星座点符号
%       Constellation_points：对应调制制式下星座点的集合
% 输出：x_ml：最大似然估计的符号
H = 1/sqrt(2)*(randn(8, 1) + 1j*randn(8, 1));     % 生成信道增益
Z = sqrt(Nt/2/(10^(snr/10)))*(randn(8, 1) + 1j*randn(8, 1));  % 模拟信道噪声
G = [x(1), x(2), x(3), 0, x(4), x(5), x(6), 0;
    -x(2)', x(1)', 0, -x(3), x(5)', -x(4)', 0, x(6);
    x(3)', 0, -x(1)', -x(2), -x(6)', 0, x(4)', x(5);
    0, -x(3)', x(2)', -x(1), 0, x(6)', -x(5)', x(4);
    -x(4), -x(5), -x(6), 0, x(1), x(2), x(3), 0;
    -x(5)', x(4)', 0, x(6), -x(2)', x(1)', 0, x(3);
    x(6)', 0, -x(4)', x(5), x(3)', 0, -x(1)', x(2);
    0, x(6)', -x(5)', -x(4), 0, x(3)', -x(2)', -x(1)];  % 生成信道发射矩阵
y = G*H + Z;    % 模拟接收信号
% 计算一些最大似然表达式中的参数
c = sum(abs(H).^2);
d = imag(H(1)*H(5)' + H(2)*H(6)' + H(3)*H(7)' - H(4)*H(8)');
f14_min = inf; f25_min = inf; f36_min = inf;
for i = 1:2^type
    x1 = Constellation_points(i);
    x2 = Constellation_points(i);
    x3 = Constellation_points(i);
    for j = 1:2^type
        x4 = Constellation_points(j);
        f14 = c*(abs(x1)^2 + abs(x4)^2) - 4*d*imag(x1*x4')...
            -2*real(x1*(H(1)*y(1)'+H(2)'*y(2)-H(3)'*y(3)-H(4)*y(4)'+H(5)*y(5)'+H(6)'*y(6)-H(7)'*y(7)-H(8)*y(8)'))...
            -2*real(x4*(-H(1)*y(5)'+H(2)'*y(6)-H(3)'*y(7)-H(4)*y(8)'+H(5)*y(1)'-H(6)'*y(2)+H(7)'*y(3)+H(8)*y(4)'));
        
        x5 = Constellation_points(j);
        f25 = c*(abs(x2)^2 + abs(x5)^2) - 4*d*imag(x2*x5')...
            -2*real(x2*(-H(1)'*y(2)+H(2)*y(1)'+H(3)'*y(4)-H(4)*y(3)'-H(5)'*y(6)+H(6)*y(5)'-H(7)'*y(8)+H(8)*y(7)'))...
            -2*real(x5*(-H(1)'*y(6)-H(2)*y(5)'-H(3)'*y(8)+H(4)*y(7)'+H(5)'*y(2)+H(6)*y(1)'-H(7)'*y(4)+H(8)*y(3)'));
        
        x6 = Constellation_points(j);
        f36 = c*(abs(x3)^2 + abs(x6)^2) - 4*d*imag(x3*x6')...
            -2*real(x3*(H(1)'*y(3)-H(2)'*y(4)+H(3)*y(1)'-H(4)*y(2)'+H(5)'*y(7)+H(6)'*y(8)+H(7)*y(5)'+H(8)*y(6)'))...
            -2*real(x6*(H(1)'*y(7)+H(2)'*y(8)-H(3)*y(5)'+H(4)*y(6)'-H(5)'*y(3)+H(6)'*y(4)+H(7)*y(1)'+H(8)*y(2)'));
        
        if f14 < f14_min
            f14_min = f14; x1_ml = x1; x4_ml = x4;
        end
        if f25 < f25_min
            f25_min = f25; x2_ml = x2; x5_ml = x5;
        end
        if f36 < f36_min
            f36_min = f36; x3_ml = x3; x6_ml = x6;
        end
    end
end
x_ml = [x1_ml, x2_ml, x3_ml, x4_ml, x5_ml, x6_ml];  % 译码信号
end