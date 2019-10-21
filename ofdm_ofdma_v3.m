%close all;
clear all;
mod_order = 16; % Oder da modulação QAM
N = 512; % Numero de sinais a ser modulado/tamaho do ODFDM symbol
%EbNo=14:2:18; % relação sinal ruido
EbNo=40:1:40;
%SNR=18;

for i=1:length(EbNo)
    
        sum_ber_u1 = 0;
        sum_ber_u2 = 0;
        sum_ber_u3 = 0;
        sum_ber_u4 = 0;
        ber_avg1_u1(i) = 0;
        ber_avg2_u2(i) = 0;
        ber_avg3_u3(i) = 0;
        ber_avg3_u4(i) = 0;

        iteracao = 0;
        SNR(i) = EbNo(i)+ 10*log10( log2(mod_order)) ;
        while(sum_ber_u1 < 200 && sum_ber_u2 < 1000)
            i
            iteracao = iteracao + 1;
            %tx 
            data1_tx(1:128) = randi([0,mod_order-1],N/4,1);
            data2_tx(1:128) = randi([0,mod_order-1],N/4,1);
            data3_tx(1:128) = randi([0,mod_order-1],N/4,1);
            data4_tx(1:128) = randi([0,mod_order-1],N/4,1);
            
            % 16QAM modulation  
            qam_symb_tx_usr1=qammod(data1_tx,mod_order, 'UnitAveragePower', true); 
            qam_symb_tx_usr2=qammod(data2_tx,mod_order, 'UnitAveragePower', true); 
            qam_symb_tx_usr3=qammod(data3_tx,mod_order, 'UnitAveragePower', true); 
            qam_symb_tx_usr4=qammod(data4_tx,mod_order, 'UnitAveragePower', true); 

            % Concatenate all users
            ofdm_symb = horzcat(qam_symb_tx_usr1,qam_symb_tx_usr2,qam_symb_tx_usr3, qam_symb_tx_usr4);
            
            %%%%%%%% Rayleigh channel(t) with multi tap %%%%%%%%
            h1 = ((randn(1,1)+1i*randn(1,1)))/ (sqrt(2) * sqrt(1));
            h2 = ((randn(1,10)+1i*randn(1,10)))/ (sqrt(2) * sqrt(10));
            h3 = ((randn(1,20)+1i*randn(1,20)))/ (sqrt(2) * sqrt(20));
            h4 = ((randn(1,30)+1i*randn(1,30)))/ (sqrt(2) * sqrt(30));
            
            %%%%%%%% Rayleigh Channel(f)
            H1 = fft(h1,(N));
            ch1_ofdm_symb= H1 .* ofdm_symb; 
            H2 = fft(h2,(N));
            ch2_ofdm_symb= H2 .* ofdm_symb;
            H3= fft(h3,(N));
            ch3_ofdm_symb= H3 .* ofdm_symb;
            H4= fft(h4,(N));
            ch4_ofdm_symb= H4 .* ofdm_symb;
            %figure; plot(1:N/3,abs(H3));

            %%%%%%%%% Adding AWGN  
            awgn_ch1_ofdm_symb= awgn(ch1_ofdm_symb , SNR(i));
            awgn_ch2_ofdm_symb= awgn(ch2_ofdm_symb , SNR(i));
            awgn_ch3_ofdm_symb= awgn(ch3_ofdm_symb , SNR(i));
            awgn_ch4_ofdm_symb= awgn(ch4_ofdm_symb , SNR(i));
            
            %rx
            data1_rx = awgn_ch1_ofdm_symb(1,1:128);
            data2_rx = awgn_ch2_ofdm_symb(1,129:256);
            data3_rx = awgn_ch3_ofdm_symb(1,257:384);
            data4_rx = awgn_ch4_ofdm_symb(1,385:512);
            
            data1_rx_demod=qamdemod((data1_rx./H1(1:128)),mod_order, 'UnitAveragePower', true);
            data2_rx_demod=qamdemod((data2_rx./H2(129:256)),mod_order, 'UnitAveragePower', true);
            data3_rx_demod=qamdemod((data3_rx./H3(257:384)),mod_order, 'UnitAveragePower', true);
            data4_rx_demod=qamdemod((data4_rx./H4(385:512)),mod_order, 'UnitAveragePower', true);

            %SER
            [SnumErrors,ser2] = symerr(data1_tx,data1_rx_demod);
            %BER USER 1
            [BnumErrors, ber1] = biterr(data1_tx,data1_rx_demod);
            %BER USER 2
            [BnumErrors, ber2] = biterr(data2_tx,data2_rx_demod);
            %BER USER 3
            [BnumErrors, ber3] = biterr(data3_tx,data3_rx_demod);
            %BER USER 4
            [BnumErrors, ber4] = biterr(data4_tx,data4_rx_demod);

            sum_ber_u1 = sum_ber_u1 + ber1;
            sum_ber_u2 = sum_ber_u2 + ber2;
            sum_ber_u3 = sum_ber_u3 + ber3;
            sum_ber_u4 = sum_ber_u4 + ber4;
            
        end 
        ber_avg1_u1(i) = sum_ber_u1/iteracao;
        ber_avg1_u2(i) = sum_ber_u2/iteracao;
        ber_avg1_u3(i) = sum_ber_u3/iteracao;
        ber_avg1_u4(i) = sum_ber_u4/iteracao;
        
        
end
%var_ebno = SNR;
%plots
% scatterplot(qam_symb_tx_f)
% grid on;
% 
% scatterplot(qam_sim_rx_f)
% grid on;

% h_test = ((randn(1,2)+1i*randn(1,2)))/ (sqrt(2) * sqrt(2));
% figure;
% plot(1:N/3,abs(fft(h_test,N/3)));


%figure('Name','Simulation Plot Window','NumberTitle','on')
figure;
semilogy(EbNo,ber_avg1_u1, '-.ob');
xlabel('EB/N0 [dB]');
ylabel('BER');
%axis([SNR(1) SNR(length(SNR)) 1e-4 .9]);
title('BER curve for 4 users using 16QAM modulation');
grid on;

hold on;
plot(EbNo,ber_avg1_u2);

hold on;
plot(EbNo,ber_avg1_u3);

hold on;
plot(EbNo,ber_avg1_u4);

legend ('BER USER1','BER USER2','BER USER3','BER USER4');

hold on;
semilogy(ebno0,ber0, '-.xk');

% figure;
% plot([1:N], ofdm_symb_t,'b.-')

                  


