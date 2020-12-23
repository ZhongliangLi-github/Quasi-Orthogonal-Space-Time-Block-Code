function x_ml = J8_coding2(x, Constellation_points, Nt, snr, type)
% 输入：x：发射的星座点符号
%       Constellation_points：对应调制制式下星座点的集合
% 输出：x_ml：最大似然估计的符号
H = 1/sqrt(2)*(randn(8, 2) + 1j*randn(8, 2));     % 生成信道增益
Z = sqrt(Nt/2/(10^(snr/10)))*(randn(8, 2) + 1j*randn(8, 2));  % 模拟信道噪声
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
d(1) = imag(H(1,1)*H(5,1)' + H(2,1)*H(6,1)' + H(3,1)*H(7,1)' - H(4,1)*H(8,1)');
d(2) = imag(H(1,2)*H(5,2)' + H(2,2)*H(6,2)' + H(3,2)*H(7,2)' - H(4,2)*H(8,2)');
f14_min = inf; f25_min = inf; f36_min = inf;
for i = 1:2^type
    x1 = Constellation_points(i);
    x2 = Constellation_points(i);
    x3 = Constellation_points(i);
    for j = 1:2^type
        x4 = Constellation_points(j);
        x5 = Constellation_points(j);
        x6 = Constellation_points(j);
        for k = 1:2
            f14(k) = c(k)*(abs(x1)^2 + abs(x4)^2) - 4*d(k)*imag(x1*x4')...
                -2*real(x1*(H(1,k)*y(1,k)'+H(2,k)'*y(2,k)-H(3,k)'*y(3,k)-H(4,k)*y(4,k)'+H(5,k)*y(5,k)'+H(6,k)'*y(6,k)-H(7,k)'*y(7,k)-H(8,k)*y(8,k)'))...
                -2*real(x4*(-H(1,k)*y(5,k)'+H(2,k)'*y(6,k)-H(3,k)'*y(7,k)-H(4,k)*y(8,k)'+H(5,k)*y(1,k)'-H(6,k)'*y(2,k)+H(7,k)'*y(3,k)+H(8,k)*y(4,k)'));
            
            f25(k) = c(k)*(abs(x2)^2 + abs(x5)^2) - 4*d(k)*imag(x2*x5')...
                -2*real(x2*(-H(1,k)'*y(2,k)+H(2,k)*y(1,k)'+H(3,k)'*y(4,k)-H(4,k)*y(3,k)'-H(5,k)'*y(6,k)+H(6,k)*y(5,k)'-H(7,k)'*y(8,k)+H(8,k)*y(7,k)'))...
                -2*real(x5*(-H(1,k)'*y(6,k)-H(2,k)*y(5,k)'-H(3,k)'*y(8,k)+H(4,k)*y(7,k)'+H(5,k)'*y(2,k)+H(6,k)*y(1,k)'-H(7,k)'*y(4,k)+H(8,k)*y(3,k)'));
            
            f36(k) = c(k)*(abs(x3)^2 + abs(x6)^2) - 4*d(k)*imag(x3*x6')...
                -2*real(x3*(H(1,k)'*y(3,k)-H(2,k)'*y(4,k)+H(3,k)*y(1,k)'-H(4,k)*y(2,k)'+H(5,k)'*y(7,k)+H(6,k)'*y(8,k)+H(7,k)*y(5,k)'+H(8,k)*y(6,k)'))...
                -2*real(x6*(H(1,k)'*y(7,k)+H(2,k)'*y(8,k)-H(3,k)*y(5,k)'+H(4,k)*y(6,k)'-H(5,k)'*y(3,k)+H(6,k)'*y(4,k)+H(7,k)*y(1,k)'+H(8,k)*y(2,k)'));
        end
        if sum(f14) < f14_min
            f14_min = sum(f14); x1_ml = x1; x4_ml = x4;
        end
        if sum(f25) < f25_min
            f25_min = sum(f25); x2_ml = x2; x5_ml = x5;
        end
        if sum(f36) < f36_min
            f36_min = sum(f36); x3_ml = x3; x6_ml = x6;
        end
    end
end
x_ml = [x1_ml, x2_ml, x3_ml, x4_ml, x5_ml, x6_ml];  % 译码信号
end