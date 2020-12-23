clc, clear, close all
%% 初始化
Nt = 8;     % 发射天线数目
Nr = 1;     % 接收天线数目
type = [1, 2, 3, 4];   % 调制方式：BPSK、QPSK、8PSK、16QAM
frame = 2000;   % 发送2000个帧
error_rate_BPSK = zeros(1,16);
error_rate_QPSK = zeros(1,16);
error_rate_8PSK = zeros(1,16);
error_rate_16QAM = zeros(1,16);
%% 循环
for snr = 0:2:30
    for k = 1:frame
        % 生成二进制比特流
        tx_bits_BPSK = randi([0,1], 1, 6*type(1));
        tx_bits_QPSK = randi([0,1], 1, 6*type(2));
        tx_bits_8PSK = randi([0,1], 1, 6*type(3));
        tx_bits_16QAM = randi([0,1], 1, 6*type(4));
        % 把二进制符号转换为四进制符号
        tx_BPSK = tx_bits_BPSK;
        tx_QPSK = 2*tx_bits_QPSK(1:2:end) + tx_bits_QPSK(2:2:end);
        tx_8PSK = 4*tx_bits_8PSK(1:3:end) + 2*tx_bits_8PSK(2:3:end) + tx_bits_8PSK(3:3:end);
        tx_16QAM = 8*tx_bits_16QAM(1:4:end) + 4*tx_bits_16QAM(2:4:end) + 2*tx_bits_16QAM(3:4:end) + tx_bits_16QAM(4:4:end);
        % 采用格雷码进行调制
        x_BPSK = pskmod(tx_BPSK, 2, 0, 'gray');     % BPSK
        x_QPSK = pskmod(tx_QPSK, 4, pi/4, 'gray');    % pi/4-QPSK
        x_8PSK = pskmod(tx_8PSK, 8, 0, 'gray');     % 8PSK
        x_16QAM = qammod(tx_16QAM, 16, 'gray', 'UnitAveragePower', true);   % 16QAM（单位化功率）
        % 不同调制制式下的星座点
        Constellation_points_BPSK = [1, -1+0*1j];
        Constellation_points_QPSK = [1+1j, -1+1j, 1-1j, -1-1j]/sqrt(2);
        Constellation_points_8PSK = [1, (1+1j)/sqrt(2), 1j, (-1+1j)/sqrt(2), -1, (-1-1j)/sqrt(2), -1j, (1-1j)/sqrt(2)];
        Constellation_points_16QAM = [1+1j, 1+3*1j, 3+1j, 3+3*1j, -1+1j, -1+3*1j, -3+1j, -3+3*1j,...
            -1-1j, -1-3*1j, -3-1j, -3-3*1j, 1-1j, 1-3*1j, 3-1j, 3-3*1j]/sqrt(10);
        % 最大似然译码
        x_ml_BPSK = J8_coding2(x_BPSK, Constellation_points_BPSK, Nt, snr, type(1));
        x_ml_QPSK = J8_coding2(x_QPSK, Constellation_points_QPSK, Nt, snr, type(2));
        x_ml_8PSK = J8_coding2(x_8PSK, Constellation_points_8PSK, Nt, snr, type(3));
        x_ml_16QAM = J8_coding2(x_16QAM, Constellation_points_16QAM, Nt, snr, type(4));
        % 对应采用格雷码解调信号
        rx_BPSK = pskdemod(x_ml_BPSK, 2, 0, 'gray');
        rx_QPSK = pskdemod(x_ml_QPSK, 4, pi/4, 'gray');
        rx_8PSK = pskdemod(x_ml_8PSK, 8, 0, 'gray');
        rx_16QAM = qamdemod(x_ml_16QAM, 16, 'gray', 'UnitAveragePower', true);
        % 转换为二进制数据
        rx_bits_BPSK = reshape(double(dec2base(rx_BPSK, 2, type(1)))'-48, 1, 6*type(1));
        rx_bits_QPSK = reshape(double(dec2base(rx_QPSK, 2, type(2)))'-48, 1, 6*type(2));
        rx_bits_8PSK = reshape(double(dec2base(rx_8PSK, 2, type(3)))'-48, 1, 6*type(3));
        rx_bits_16QAM = reshape(double(dec2base(rx_16QAM, 2, type(4)))'-48, 1, 6*type(4));
        % 统计错误比特数
        error_sum_BPSK(k) = sum(abs(rx_bits_BPSK - tx_bits_BPSK));
        error_sum_QPSK(k) = sum(abs(rx_bits_QPSK - tx_bits_QPSK));
        error_sum_8PSK(k) = sum(abs(rx_bits_8PSK - tx_bits_8PSK));
        error_sum_16QAM(k) = sum(abs(rx_bits_16QAM - tx_bits_16QAM));
    end
    error_rate_BPSK(snr/2+1) = sum(error_sum_BPSK)/(6*type(1)*frame);
    error_rate_QPSK(snr/2+1) = sum(error_sum_QPSK)/(6*type(2)*frame);
    error_rate_8PSK(snr/2+1) = sum(error_sum_8PSK)/(6*type(3)*frame);
    error_rate_16QAM(snr/2+1) = sum(error_sum_16QAM)/(6*type(4)*frame);
end
%% 作图
snr = 0:2:30;
semilogy(snr, error_rate_BPSK, 'LineWidth', 2)
hold on
semilogy(snr, error_rate_QPSK, 'LineWidth', 2)
semilogy(snr, error_rate_8PSK, 'LineWidth', 2)
semilogy(snr, error_rate_16QAM, 'LineWidth', 2)
xlabel('信噪比[dB]');ylabel('误码率BER');title('信噪比和误码率的关系')
legend('BPSK', 'QPSK', '8PSK', '16QAM')
grid on