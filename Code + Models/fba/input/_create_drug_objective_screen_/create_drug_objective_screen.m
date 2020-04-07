clear, clc, close all

% drug to choose
drug = 'dox';
drug_choices = {'5fu','cis','cpa','dox'};

% list of reactions
rxns = {
    {'5fu_transport_e_c','5fu_transport_c_e','5fu_fdurd','fdurd_fdump','fdudp_fdump','5fu_furd','furd_fump','5fu_fump','fump_fudp','fudp_fdudp','5fu_dhfu','dhfu_fupa','fupa_fbal','fbal_transport_c_e'},
    {'cis_transport_e_c','cis_transport_c_e','cisGSH_transport_c_e','cisMT_transport_c_e','cis_transport_c_m','cis_cis-mtDNA','cis_cisGSH','cis_cisMT'},
    {'cpa_transport_e_c','cpa_transport_c_e','dce-cpa_transport_c_e','caa_transport_c_e','4k-cpa_transport_c_e','alcph_transport_c_e','cph_transport_c_e','phmust-GSH_transport_c_e','acr-GSH_transport_c_e','cpa_4oh-cpa','4oh-cpa_aldph','aldph_phmust','cpa_dce-cpa','4oh-cpa_4k-cpa','aldph_alcph','aldph_cph','phmust_phmust-GSH','acr_acr-GSH'},
    {'doxQ_transport_e_c','doxQ_transport_c_e_noATP','doxQ_transport_c_e_withATP','doxHQ_transport_c_e','doxdeoxy_transport_c_e','daun_transport_c_e','doxQ_doxSQ_NOS3','doxQ_doxSQ_POR','doxSQ_doxQ','doxQ_doxHQ','doxQ_doxdeoxy'}
 };

% list of metabolites
mets = {
    {'5fu[c]','fdurd[c]','fdump[c]','fdudp[c]','furd[c]','fump[c]','fudp[c]','dhfu[c]','fupa[c]','fbal[c]','atp[c]','adp[c]','2dr1p[c]','r1p[c]','prpp[c]','trdrd[c]','trdox[c]','nadph[c]','nadp[c]'},
    {'cis[c]','cis[m]','cis-mtDNA[m]','cisGSH[c]','cisMT[c]','atp[c]','adp[c]','o2s[m]','gthrd[c]','gthox[c]'},
    {'cpa[c]','dce-cpa[c]','caa[c]','4oh-cpa[c]','4k-cpa[c]','aldph[c]','alcph[c]','cph[c]','phmust[c]','acr[c]','phmust-GSH[c]','acr-GSH[c]','atp[c]','adp[c]','nadph[c]','nadp[c]','nad[c]','nadh[c]','gthrd[c]'},
    {'doxQ[c]','doxSQ[c]','doxHQ[c]','doxdeoxy[c]','daun[c]','atp[c]','adp[c]','arg_L[c]','nadph[c]','citr_L[c]','nadp[c]','o2s[c]'}
};

% list of other reactions
others = {
    {'1 nadph[c] --> 1 nadp[c] + 1 h[c]','1 nadp[c] + 1 h[c] --> 1 nadph[c]','1 atp[c] + 1 h2o[c] --> 1 adp[c] + 1 pi[c] + 1 h[c]','1 adp[c] + 1 pi[c] + 1 h[c] --> 1 atp[c] + 1 h2o[c]','1 nadp[c] + 1 trdrd[c] --> 1 h[c] + 1 nadph[c] + 1 trdox[c]','1 h[c] + 1 nadph[c] + 1 trdox[c] --> 1 nadp[c] + 1 trdrd[c]'},
    {'1 atp[c] + 1 h2o[c] --> 1 adp[c] + 1 pi[c] + 1 h[c]','1 adp[c] + 1 pi[c] + 1 h[c] --> 1 atp[c] + 1 h2o[c]','1 o2[m] --> 1 o2s[m]','1 o2s[m] --> 1 o2[m]','1 gthox[c] + 1 h[c] + 1 nadph[c] --> 2 gthrd[c] + 1 nadp[c]','2 gthrd[c] + 1 nadp[c] --> 1 gthox[c] + 1 h[c] + 1 nadph[c]'},
    {'1 atp[c] + 1 h2o[c] --> 1 adp[c] + 1 pi[c] + 1 h[c]','1 adp[c] + 1 pi[c] + 1 h[c] --> 1 atp[c] + 1 h2o[c]','1 nadph[c] --> 1 nadp[c] + 1 h[c]','1 nadp[c] + 1 h[c] --> 1 nadph[c]','1 nadh[c] --> 1 nad[c] + 1 h[c]','1 nad[c] + 1 h[c] --> 1 nadh[c]','1 gthox[c] + 1 h[c] + 1 nadph[c] --> 2 gthrd[c] + 1 nadp[c]','2 gthrd[c] + 1 nadp[c] --> 1 gthox[c] + 1 h[c] + 1 nadph[c]'},
    {'1 atp[c] + 1 h2o[c] --> 1 adp[c] + 1 pi[c] + 1 h[c]','1 adp[c] + 1 pi[c] + 1 h[c] --> 1 atp[c] + 1 h2o[c]','1 nadph[c] --> 1 nadp[c] + 1 h[c]','1 nadp[c] + 1 h[c] --> 1 nadph[c]','1 o2[m] --> 1 o2s[m]','1 o2s[m] --> 1 o2[m]'}
};

% only run if input folder doesn't already exist
if ~exist(sprintf('../%s/',drug))

	% create output folder
	mkdir(sprintf('../%s/',drug));
	
	% load Recon3
	load('../../../data/recon/recon3d_qflux.mat')

	% iterate over reactions
	for i = 1:length(rxns{strcmp(drug_choices,drug)})

        % make copy of original file
        copyfile(sprintf('%s.xlsx',drug),sprintf('../%s/rxn_%d.xlsx',drug,i));

        % edit file
        xlwrite(sprintf('../%s/rxn_%d.xlsx',drug,i),{'X'},'Objective Function','A2');
        xlwrite(sprintf('../%s/rxn_%d.xlsx',drug,i),{'MAX'},'Objective Function','B2');
        xlwrite(sprintf('../%s/rxn_%d.xlsx',drug,i),{'1'},'Objective Function','C2');
        xlwrite(sprintf('../%s/rxn_%d.xlsx',drug,i),{rxns{strcmp(drug_choices,drug)}{i}},'Objective Function','D2');
    end
    
    % iterate over metabolites
	for i = 1:length(mets{strcmp(drug_choices,drug)})

        % make copy of original file
        copyfile(sprintf('%s.xlsx',drug),sprintf('../%s/met_%d.xlsx',drug,i));

        % edit file
        xlwrite(sprintf('../%s/met_%d.xlsx',drug,i),{'X'},'Objective Function','A2');
        xlwrite(sprintf('../%s/met_%d.xlsx',drug,i),{'MAX'},'Objective Function','B2');
        xlwrite(sprintf('../%s/met_%d.xlsx',drug,i),{'1'},'Objective Function','C2');
        xlwrite(sprintf('../%s/met_%d.xlsx',drug,i),{sprintf('1 %s --> ',mets{strcmp(drug_choices,drug)}{i})},'Objective Function','D2');
    end
    
    % iterate over other reactions
	for i = 1:length(others{strcmp(drug_choices,drug)})

        % make copy of original file
        copyfile(sprintf('%s.xlsx',drug),sprintf('../%s/other_%d.xlsx',drug,i));

        % edit file
        xlwrite(sprintf('../%s/other_%d.xlsx',drug,i),{'X'},'Objective Function','A2');
        xlwrite(sprintf('../%s/other_%d.xlsx',drug,i),{'MAX'},'Objective Function','B2');
        xlwrite(sprintf('../%s/other_%d.xlsx',drug,i),{'1'},'Objective Function','C2');
        xlwrite(sprintf('../%s/other_%d.xlsx',drug,i),{others{strcmp(drug_choices,drug)}{i}},'Objective Function','D2');
	end

% raise error if input folder already exists
else
	error('Input folder already exists');
end