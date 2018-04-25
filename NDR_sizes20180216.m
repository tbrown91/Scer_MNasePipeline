load ../gene_annotations/harry_annotations.mat;
load ../gene_annotations/RPG_geneCoords.mat;
load ../gene_annotations/OX_geneCoords.mat;
load ../gene_annotations/RB_geneCoords.mat;
load ../gene_annotations/RC_geneCoords.mat;
load ../gene_annotations/harry_5Clusters_1.mat;
load ../gene_annotations/harry_5Clusters_2.mat;
load ../gene_annotations/harry_5Clusters_3.mat;
load ../gene_annotations/harry_5Clusters_4.mat;
load ../gene_annotations/harry_5Clusters_5.mat;

load ../MNase_files/DANPOS_nucPositions_BY_D.mat;
load ../MNase_files/DANPOS_nucPositions_BY_15mG.mat;
load ../MNase_files/DANPOS_nucPositions_BY_60mG.mat;

gene_datasets = cell(10,1);
gene_datasets{1}=harry_annotations;
gene_datasets{2}=RPG_geneCoords;
gene_datasets{3}=OX_geneCoords;
gene_datasets{4}=RB_geneCoords;
gene_datasets{5}=RC_geneCoords;
gene_datasets{6}=harry_5Clusters_1;
gene_datasets{7}=harry_5Clusters_2;
gene_datasets{8}=harry_5Clusters_3;
gene_datasets{9}=harry_5Clusters_4;
gene_datasets{10}=harry_5Clusters_5;

dataset_labels = {'All genes','RPGs','OX genes','RB genes','RC genes','Harry cluster 1','Harry cluster 2','Harry cluster 3','Harry cluster 4','Harry cluster 5'};

clear harry_annotations;
clear RPG_geneCoords;
clear OX_geneCoords;
clear RB_geneCoords;
clear RC_geneCoords;
clear harry_5Clusters_1;
clear harry_5Clusters_2;
clear harry_5Clusters_3;
clear harry_5Clusters_4;
clear harry_5Clusters_5;
%% Find -1 and +1 nucleosomes, report size of NDR

TSS_start=-20;

nuc_positions=cell(length(gene_datasets),1);

num_genes=zeros(length(gene_datasets),1);
for d=1:length(num_genes)
    for chr=1:size(gene_datasets{d},1)
        for strand=1:2
            num_genes(d)=num_genes(d)+size(gene_datasets{d}{chr,strand},1);
        end
    end
end
NDR_sizes = cell(length(gene_datasets),1);

for d=1:length(gene_datasets)
    NDR_sizes{d}=nan(3,num_genes(d));
    nuc_positions{d} = nan(3,2,num_genes(d));
    k=1;
    for chr=1:16
        strand=1;
        for gene=1:size(gene_datasets{d}{chr,strand},1)
            TSS=gene_datasets{d}{chr,strand}(gene,1);
            %Find -1,+1,+2,+3 nucleosomes relative to each TSS for each
            %condition
            idx_D = find(DANPOS_nucPositions_BY_D{chr,strand}(TSS-500:TSS+1000))-500;
            idx_15mG = find(DANPOS_nucPositions_BY_15mG{chr,strand}(TSS-500:TSS+1000))-500;
            idx_60mG = find(DANPOS_nucPositions_BY_60mG{chr,strand}(TSS-500:TSS+1000))-500;
            
            try
                nuc_positions{d}(1,1,k)=max(idx_D(idx_D<=TSS_start-50));
            end
            try
                nuc_positions{d}(2,1,k)=max(idx_15mG(idx_15mG<=TSS_start-50));
            end
            try
                nuc_positions{d}(3,1,k)=max(idx_60mG(idx_60mG<=TSS_start-50));
            end
            
            try
                nuc_positions{d}(1,2,k)=min(idx_D((idx_D>=TSS_start)&(idx_D<=TSS_start+100)));
            end
         
            try
                nuc_positions{d}(2,2,k)=min(idx_15mG((idx_15mG>=TSS_start)&(idx_15mG<=TSS_start+100)));
            end
          
            
            try
                nuc_positions{d}(3,2,k)=min(idx_60mG((idx_60mG>=TSS_start)&(idx_60mG<=TSS_start+100)));
            end          

            k=k+1;
        end
        strand=2;
        for gene=1:size(gene_datasets{d}{chr,strand},1)
            TSS=gene_datasets{d}{chr,strand}(gene,1);
            %Find -1,+1,+2,+3 nucleosomes relative to each TSS for each
            %condition
            idx_D = fliplr(-(find(DANPOS_nucPositions_BY_D{chr,strand}(TSS-1000:TSS+500))-1000)');
            idx_15mG = fliplr(-(find(DANPOS_nucPositions_BY_15mG{chr,strand}(TSS-1000:TSS+500))-1000)');
            idx_60mG = fliplr(-(find(DANPOS_nucPositions_BY_60mG{chr,strand}(TSS-1000:TSS+500))-1000)');
            
            try
                nuc_positions{d}(1,1,k)=max(idx_D(idx_D<=TSS_start-50));
            end
            try
                nuc_positions{d}(2,1,k)=max(idx_15mG(idx_15mG<=TSS_start-50));
            end
            try
                nuc_positions{d}(3,1,k)=max(idx_60mG(idx_60mG<=TSS_start-50));
            end
            
            try
                nuc_positions{d}(1,2,k)=min(idx_D((idx_D>=TSS_start)&(idx_D<=TSS_start+100)));
            end
                       
            try
                nuc_positions{d}(2,2,k)=min(idx_15mG((idx_15mG>=TSS_start)&(idx_15mG<=TSS_start+100)));
            end
            
            try
                nuc_positions{d}(3,2,k)=min(idx_60mG((idx_60mG>=TSS_start)&(idx_60mG<=TSS_start+100)));
            end
            
            k=k+1;
        end
    end
    p_vals = nan(2,1);
    for t=1:3
        NDR_sizes{d}(t,:)=reshape(nuc_positions{d}(t,2,:)-nuc_positions{d}(t,1,:),1,num_genes(d));
        
    end
    [~,p_vals(1)] = ttest(NDR_sizes{d}(2,:)-NDR_sizes{d}(1,:));
    [~,p_vals(2)] = ttest(NDR_sizes{d}(3,:)-NDR_sizes{d}(1,:));
    
    
    fig=figure('position',[100 100 1080 800]);
    boxplot(NDR_sizes{d}','orientation','horizontal');
    text(200,2.3,['p = ',num2str(p_vals(1))],'fontsize',16);
    text(200,3.3,['p = ',num2str(p_vals(2))],'fontsize',16);
    xlim([0 500]);
    xlabel('NDR Size (bp)','fontsize',24);
    ax=gca;
    ax.YTickLabel = {'D','15mG','60mG'};
    ax.FontSize=24;
    stitle = suptitle(dataset_labels{d});
    stitle.FontSize=24;
    %saveas(fig,strcat('./figures_20180216/',dataset_labels{d},'_NDRSizes.fig'));
end
