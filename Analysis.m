clearvars;
close all;

cell_id = 135:269; % all of the cell ID numbers
isi_avg = []; %to put all of the ISI averages per file


for input_file_number = 2:12 %Need to run code for each file. This for loop runs the file for 10 times
number_str = num2str(input_file_number);
spike_file = readmatrix(strcat('Gamma',number_str,'.txt')); %Loads file (first column is times and second column is ID number

%put code here

isis_spk = []; % column vector for storing the ISI value for each cell to be averaged later
isi_rand = [];
rand_id = randi([135 269],1,20);
sample_isi_avg = [];


% find spikes corresponding to each pyramidal cell
cell_spikes = cell(length(cell_id),1); % matrix of cells where each item contains the list of spikes for a given cell
for idx=1:length(cell_id)
    current_id = cell_id(idx); % the current ID
    sel = spike_file(:,2)==current_id; % select indicies within spike_file that correspond to the current ID
    current_times = spike_file(sel,1); % select the times of all spikes that correspond to the current ID
    curr_isi = diff(current_times); % calculate isi values
    
    isis_spk = [isis_spk;curr_isi]; %appends to isis_spk the recently calculated isi
    if any(rand_id(:) == current_id) 
        isi_rand = [isi_rand, curr_isi'];
    end
    
    
end

isi_avg = [isi_avg, mean(isis_spk)]; 
sample_isi_avg = [sample_isi_avg, mean(isi_rand)];


end

disp(isi_avg);
disp(sample_isi_avg);