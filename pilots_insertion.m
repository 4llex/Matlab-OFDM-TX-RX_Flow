function [ ofdm_symbol_with_pilots ] = pilots_insertion(Nutil, ofdm_symb, Pilots)
N = Nutil;% 448;
%ofdm_symb = zeros(1,N);
P = Pilots;%64;

f = 1;
ofdm_symbol_with_pilots= [];
    j = 1;
    k = 1;
    g = 1;
    np=1;
    while g < (N/7)+P+1
        if j==1
            ofdm_symbol_with_pilots(f) = 1;
            j = 2;
            f = f+1;
            np = np+1;
        else
           ofdm_symbol_with_pilots(f:f+6) = ofdm_symb(k:k+6);
           j = 1;
           k = k+7;
           f = f+7;
        end
        g = g+1;
    end