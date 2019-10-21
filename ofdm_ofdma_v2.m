%close all;
clear all;
mod_order = 16; % Oder da modulação QAM
N = 300; % Numero de sinais a ser modulado/tamaho do ODFDM symbol
%EbNo=14:2:18; % relação sinal ruido
EbNo=1:2:6;
%SNR=18;

for i=1:length(EbNo)
    
        soma_ber_user1 = 0;
        soma_ber_user2 = 0;
        soma_ber_user3 = 0;
        var_ber_media1(i) = 0;
        var_ber_media2(i) = 0;
        var_ber_media3(i) = 0;

        iteracao = 0;
        SNR(i) = EbNo(i)+ 10*log10( log2(mod_order)) ;
        while(soma_ber_user1 < 400)
            i
            iteracao = iteracao + 1;
            %tx 
            data1_tx(1:100) = randi([0,mod_order-1],N/3,1);
            data2_tx(1:100) = randi([0,mod_order-1],N/3,1);
            data3_tx(1:100) = randi([0,mod_order-1],N/3,1);
            
            % 16QAM modulation  
            qam_symb_tx_usr1=qammod(data1_tx,mod_order, 'UnitAveragePower', true); 
            qam_symb_tx_usr2=qammod(data2_tx,mod_order, 'UnitAveragePower', true); 
            qam_symb_tx_usr3=qammod(data3_tx,mod_order, 'UnitAveragePower', true); 

            
            %%%%%%%% Rayleigh channel with multi tap %%%%%%%%
            h1 = ((randn(1,1)+1i*randn(1,1)))/ (sqrt(2) * sqrt(1));
            h2 = ((randn(1,10)+1i*randn(1,10)))/ (sqrt(2) * sqrt(10));
            h3 = ((randn(1,20)+1i*randn(1,20)))/ (sqrt(2) * sqrt(20));
            
            H1 = fft(h1,(N/3));
            %figure; plot(1:N/3,abs(H1));
            rayleigh_symb_usr1= H1 .* qam_symb_tx_usr1; 
            
            H2 = fft(h2,(N/3));
            %figure; plot(1:N/3,abs(H2));
            rayleigh_symb_usr2= H2 .* qam_symb_tx_usr2;
            
            H3= fft(h3,(N/3));
            %figure; plot(1:N/3,abs(H3));
            rayleigh_symb_usr3= H3 .* qam_symb_tx_usr3;

            
            %Adding AWGN / tem que ser na freq?
            %ofdm_symb_awgn_t=ifft(awgn( fft(ofdm_symb_t) , SNR, 'measured'));
            ofdm_symb = horzcat(rayleigh_symb_usr1,rayleigh_symb_usr2,rayleigh_symb_usr3);
            ofdm_symb_awgn= awgn(ofdm_symb , SNR(i));%, 'measured');
            
            %rx
            data1_rx = ofdm_symb_awgn(1,1:100);
            data2_rx = ofdm_symb_awgn(1,101:200);
            data3_rx = ofdm_symb_awgn(1,201:300);
            
            data1_rx_demod=qamdemod((data1_rx./H1),mod_order, 'UnitAveragePower', true);
            data2_rx_demod=qamdemod((data2_rx./H2),mod_order, 'UnitAveragePower', true);
            data3_rx_demod=qamdemod((data3_rx./H3),mod_order, 'UnitAveragePower', true);
            
            %SER
            [SnumErrors,ser2] = symerr(data1_tx,data1_rx_demod);
            %BER USER 1
            [BnumErrors, ber1] = biterr(data1_tx,data1_rx_demod);
            %BER USER 2
            [BnumErrors, ber2] = biterr(data2_tx,data2_rx_demod);
            %BER USER 3
            [BnumErrors, ber3] = biterr(data3_tx,data3_rx_demod);


            soma_ber_user1 = soma_ber_user1 + ber1;
            soma_ber_user2 = soma_ber_user2 + ber2;
            soma_ber_user3 = soma_ber_user3 + ber3;
            
            
        end 
        var_ber_media1(i) = soma_ber_user1/iteracao;
        var_ber_media2(i) = soma_ber_user2/iteracao;
        var_ber_media3(i) = soma_ber_user3/iteracao;
        
        
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
semilogy(EbNo,var_ber_media1, '-.ob');
xlabel('EB/N0 [dB]');
ylabel('BER');
%axis([SNR(1) SNR(length(SNR)) 1e-4 .9]);
title('Bit error probability curve for 16QAM modulation');
grid on;

hold on;
plot(EbNo,var_ber_media2);

hold on;
plot(EbNo,var_ber_media3);

legend ('BER USER1','BER USER2','BER USER3');

%hold on;
%semilogy(ebno0,ber0, '-.xk');

% figure;
% plot([1:N], ofdm_symb_t,'b.-')

                  


