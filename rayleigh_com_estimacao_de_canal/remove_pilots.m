function [ data_without_pilots ] = remove_pilots( pilots_index, first_index, data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    idx = pilots_index;
    ofdm_symb_with_pilot = data;
    ft = first_index;
    
    for i=1:length(idx)
        ofdm_symb_with_pilot(idx(i)) = [];
        idx = idx - 1;
    end
    
    data_without_pilots = ofdm_symb_with_pilot;
end

