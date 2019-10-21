%close all;
clear all;
mod_order = 16; % Oder da modulação QAM
N = 512; % Numero de sinais a ser modulado/tamaho do OFDM symbol
P = 64; 
Nutil = N-P;
%EbNo=14:2:18; % relação sinal ruido
EbNo=1:2:30;
%SNR=18;
tic;



for i=1:1:length(EbNo)
    
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
        
        while(sum_ber_u1 < 300 && sum_ber_u2 < 700)
            i
            iteracao = iteracao + 1;
            %tx 
            data1_tx(1:Nutil/4) = randi([0,mod_order-1],Nutil/4,1);
            data2_tx(1:Nutil/4) = randi([0,mod_order-1],Nutil/4,1);
            data3_tx(1:Nutil/4) = randi([0,mod_order-1],Nutil/4,1);
            data4_tx(1:Nutil/4) = randi([0,mod_order-1],Nutil/4,1);
            
            % 16QAM modulation  
            qam_symb_tx_usr1=qammod(data1_tx,mod_order, 'UnitAveragePower', true); 
            qam_symb_tx_usr2=qammod(data2_tx,mod_order, 'UnitAveragePower', true); 
            qam_symb_tx_usr3=qammod(data3_tx,mod_order, 'UnitAveragePower', true); 
            qam_symb_tx_usr4=qammod(data4_tx,mod_order, 'UnitAveragePower', true); 

            % Concatenate all users
            ofdm_symb = horzcat(qam_symb_tx_usr1,qam_symb_tx_usr2,qam_symb_tx_usr3, qam_symb_tx_usr4);
            
            % Criando Símbolo Piloto
            pilot_symbol = pilots_insertion(Nutil, zeros(1, N-P), P);
            %pilot_symbol = zeros(1, N-P);
            
            % Pilot insertion 
            ofdm_symb_with_pilot= pilots_insertion(Nutil, ofdm_symb, P);
            
            %%%%%%%% Rayleigh channel(t) with multi tap %%%%%%%%
            h1 = ((randn(1,1)+1i*randn(1,1)))/ (sqrt(2) * sqrt(1));
            h2 = ((randn(1,10)+1i*randn(1,10)))/ (sqrt(2) * sqrt(10));
            h3 = ((randn(1,20)+1i*randn(1,20)))/ (sqrt(2) * sqrt(20));
            h4 = ((randn(1,30)+1i*randn(1,30)))/ (sqrt(2) * sqrt(30));
            
            %%%%%%%% Rayleigh Channel(f)
            H1 = fft(h1,(N));
            ch1_ofdm_symb= H1 .* ofdm_symb_with_pilot; 
            H2 = fft(h2,(N));
            ch2_ofdm_symb= H2 .* ofdm_symb_with_pilot;
            H3= fft(h3,(N));
            ch3_ofdm_symb= H3 .* ofdm_symb_with_pilot;
            H4= fft(h4,(N));
            ch4_ofdm_symb= H4 .* ofdm_symb_with_pilot;
            %figure; plot(1:N,abs(H4));

            %%%%%%%%% Adding AWGN  
            awgn_ch1_ofdm_symb= awgn(ch1_ofdm_symb , SNR(i));
            awgn_ch2_ofdm_symb= awgn(ch2_ofdm_symb , SNR(i));
            awgn_ch3_ofdm_symb= awgn(ch3_ofdm_symb , SNR(i));
            awgn_ch4_ofdm_symb= awgn(ch4_ofdm_symb , SNR(i));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%% RX - Receiver %%%%%%%%%%%%%%%%%%%% 
            %%% Getting pilots for each channel 
            pilotos_rx_1 = awgn_ch1_ofdm_symb .* pilot_symbol;
            pilotos_rx_2 = awgn_ch2_ofdm_symb .* pilot_symbol;
            pilotos_rx_3 = awgn_ch3_ofdm_symb .* pilot_symbol;
            pilotos_rx_4 = awgn_ch4_ofdm_symb .* pilot_symbol;
            
            %%% Get index of the pilots 
            pilot_positions = find(abs(pilotos_rx_1)); %% The positions are the same for all channels
            
            %%% Estimate the channels for each user
            Hestimado_1 = interpft((pilotos_rx_1(pilot_positions)),N);
            Hestimado_2 = interpft((pilotos_rx_2(pilot_positions)),N);
            Hestimado_3 = interpft((pilotos_rx_3(pilot_positions)),N);
            Hestimado_4 = interpft((pilotos_rx_4(pilot_positions)),N);
            
%                                         figure; plot(1:N,abs(H4), '-.');
%                                         hold on
%                                         plot(1:N,interpft(abs(pilotos_rx_4(pilot_positions)),N), '--.r');
%                                         axis([0 N 0 3]); 
            
            %%% Users receiveing data
            data1_rx = awgn_ch1_ofdm_symb;
            data2_rx = awgn_ch2_ofdm_symb;
            data3_rx = awgn_ch3_ofdm_symb;
            data4_rx = awgn_ch4_ofdm_symb;
            
            %%% Equalização de canal
            data1_rx_demod=qamdemod((data1_rx./Hestimado_1),mod_order, 'UnitAveragePower', true);
            data2_rx_demod=qamdemod((data2_rx./Hestimado_2),mod_order, 'UnitAveragePower', true);
            data3_rx_demod=qamdemod((data3_rx./Hestimado_3),mod_order, 'UnitAveragePower', true);
            data4_rx_demod=qamdemod((data4_rx./Hestimado_4),mod_order, 'UnitAveragePower', true);
            
            %%% Removing Pilots and Getting usefull data for each user
            data1_demod_useful = remove_pilots(pilot_positions, 1 , data1_rx_demod);
            data2_demod_useful = remove_pilots(pilot_positions, 1 , data2_rx_demod);
            data3_demod_useful = remove_pilots(pilot_positions, 1 , data3_rx_demod);
            data4_demod_useful = remove_pilots(pilot_positions, 1 , data4_rx_demod);

            %SER
            %[SnumErrors,ser2] = symerr(data1_tx,data1_rx_demod);
            %BER USER 1
            [BnumErrors1, ber1] = biterr(data1_tx,data1_demod_useful(1:112));
            %BER USER 2
            [BnumErrors2, ber2] = biterr(data2_tx,data2_demod_useful(113:224));
            %BER USER 3
            [BnumErrors3, ber3] = biterr(data3_tx,data3_demod_useful(225:336));
            %BER USER 4
            [BnumErrors4, ber4] = biterr(data4_tx,data4_demod_useful(337:448));

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
toc

%figure('Name','Simulation Plot Window','NumberTitle','on')
figure;
semilogy(EbNo,ber_avg1_u1, '-.xb');
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
semilogy(ebno0,ber0, '-.ok');

% figure;
% plot([1:N], ofdm_symb_t,'b.-')

                  


