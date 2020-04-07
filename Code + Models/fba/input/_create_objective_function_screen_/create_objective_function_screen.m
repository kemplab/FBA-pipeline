% user inputs
input_filename = '190308_Vandy_Melanoma_Skin.xlsx';
output_folder_name = '190308_Vandy_Melanoma_Skin';
output_file_name = 'met';

% only run if input folder doesn't already exist
if ~exist(sprintf('../%s/',output_folder_name))

	% create output folder
	mkdir(sprintf('../%s/',output_folder_name));
	
	% load Recon3
	load('../../../data/recon/recon3d_qflux.mat')
	
	% list of all metabolites
	list_of_metabolites = {};
	for i = 1:length(model.mets)
	
		% if has KEGG ID
		if ~strcmp(model.metKEGGID{i},'')
			list_of_metabolites{end+1} = model.mets{i}(1:end-3);
		end
	end
	list_of_metabolites = unique(list_of_metabolites);
	
	% iterate over unique metabolites
	for i = 1:length(list_of_metabolites)
	
		% if not 5FU metabolite
		if ~any(strcmp({'5fu','fdurd','fdump','fdudp','furd','fump','fudp','dhfu','fupa','fbal'},list_of_metabolites{i}))

			% make copy of original file
			copyfile(input_filename,sprintf('../%s/%s_%s.xlsx',output_folder_name,output_file_name,list_of_metabolites{i}));

			% edit file
			xlwrite(sprintf('../%s/%s_%s.xlsx',output_folder_name,output_file_name,list_of_metabolites{i}),{'X'},'Objective Function','A2');
			xlwrite(sprintf('../%s/%s_%s.xlsx',output_folder_name,output_file_name,list_of_metabolites{i}),{'MAX'},'Objective Function','B2');
			xlwrite(sprintf('../%s/%s_%s.xlsx',output_folder_name,output_file_name,list_of_metabolites{i}),{sprintf('all(1 %s[] --> )',list_of_metabolites{i})},'Objective Function','C2');
		end
	end

% raise error if input folder already exists
else
	error('Input folder already exists');
end


