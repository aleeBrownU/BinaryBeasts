clearvars;
close all;

cell_id = 135:269;



for input_file_number = 2:12
number_str = num2str(input_file_number);
spike_times = readmatrix(strcat('Gamma',number_str,'.txt'));



end