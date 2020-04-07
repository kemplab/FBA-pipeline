clear, clc, close all

% user inputs
input_filename = 'dox.xlsx';
output_folder_name = 'ko_dox';
output_file_name = 'gene';
knockdown_fraction = 0;

% only run if input folder doesn't already exist
if ~exist(sprintf('../%s/',output_folder_name))

	% create output folder
	mkdir(sprintf('../%s/',output_folder_name));
	
	% create excel files
	for i = 0:3268

		% make copy of file
		copyfile(input_filename,sprintf('../%s/%s_%d.xlsx',output_folder_name,output_file_name,i));

		% edit file
		if i >= 1
			xlwrite(sprintf('../%s/%s_%d.xlsx',output_folder_name,output_file_name,i),'X','Gene Knockdown',sprintf('a%d:a%d',i+2,i+2));
			xlwrite(sprintf('../%s/%s_%d.xlsx',output_folder_name,output_file_name,i),knockdown_fraction,'Gene Knockdown',sprintf('b%d:b%d',i+2,i+2));
		end
	end

% raise error if input folder already exists
else
	error('Input folder already exists');
end


