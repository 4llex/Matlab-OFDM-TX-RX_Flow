%close all;
clear all;
mod_order = 16; % Oder da modulação QAM
N = 16; % Numero de sinais a ser modulado/tamaho do ODFDM symbol
%EbNo=14:2:18; % relação sinal ruido
EbNo=1:2:21;
%SNR=18;
fading= 'yes';

for i=1:length(EbNo)
    
        var_ber_soma = 0;
        var_ber(i) = 0;
        iteracao = 0;
        SNR(i)     = EbNo(i)+ 10*log10( log2(mod_order)) ;
        while(var_ber_soma < 100)
            i
            iteracao = iteracao + 1
            %tx
            data_tx = randi([0,mod_order-1],N,1);
            qam_symb_tx_f=qammod(data_tx,mod_order, 'UnitAveragePower', true); 
            %ofdm_symb_t=ifft(qam_symb_tx_f, N);
            
            if strcmpi(fading,'yes')
                %Rayleigh channel
                h = (randn+1i*randn)/sqrt(2);
                H=fft(h, N).';   
                rayleigh_ofdm_symb= H .* qam_symb_tx_f; 
            elseif strcmpi(fading,'no')         
                h = 1;
                H=fft(h, N).';   
                rayleigh_ofdm_symb= H .* qam_symb_tx_f; 
            end
            
            %Adding AWGN / tem que ser na freq?
            %ofdm_symb_awgn_t=ifft(awgn( fft(ofdm_symb_t) , SNR, 'measured'));
            ofdm_symb_awgn= awgn(rayleigh_ofdm_symb , SNR(i));%, 'measured');
            
            %rx
            data_rx=qamdemod((ofdm_symb_awgn./H),mod_order, 'UnitAveragePower', true);
            %data_rx=qamdemod((ofdm_symb_awgn),mod_order, 'UnitAveragePower', true);
            
            %SER
            ser= sum(data_tx~=data_rx);
            [SnumErrors,ser2] = symerr(data_tx,data_rx);
            %BER
            [BnumErrors, ber] = biterr(data_tx,data_rx);

            var_ber_soma = var_ber_soma + ber
            %var_ser(i) = ser2;
        end 
        var_ber(i) = var_ber_soma/iteracao;
        
        
end
%var_ebno = SNR;
%plots
% scatterplot(qam_symb_tx_f)
% grid on;
% 
% scatterplot(qam_sim_rx_f)
% grid on;

% h_test = ((randn(1,2)+1i*randn(1,)))/ (sqrt(2) * sqrt(3));
% figure;
% plot([1:128],abs(fft(h_test,128)));
% 
% h_test = ((randn(1,1)+1i*randn(1,1)))/ (sqrt(2) * sqrt(1))
% figure;
% plot([1:128],abs(fft(h_test,128)));

%figure('Name','Simulation Plot Window','NumberTitle','on')
figure;
semilogy(EbNo,var_ber, '-.ob');
xlabel('EB/N0 [dB]');
ylabel('BER');
%axis([SNR(1) SNR(length(SNR)) 1e-4 .9]);
title('Bit error probability curve for 16QAM modulation');
legend ('BER');
grid on;


% hold on;
% semilogy(ebno0,ber0, '-.xk');


% figure;
% plot([1:N], ofdm_symb_t,'b.-')


% figure;
% [psd,f] = periodogram(ofdm_symb_t, rectwin(length(ofdm_symb_t)), N*2, ...,
%                       1);
                  
% plot(f,10*log10(psd));
% title('Power Spectral Density OFDM - Periodogram');      
% xlabel('Frequency (Hz)')         
% ylabel('Power');

