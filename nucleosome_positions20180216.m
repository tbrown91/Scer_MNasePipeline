load ./harry_annotations.mat;
load ./RPG_geneCoords.mat;
load ./OX_geneCoords.mat;
load ./RB_geneCoords.mat;
load ./RC_geneCoords.mat;

load ./DANPOS_nucPositions_BY_D.mat;
load ./DANPOS_nucPositions_BY_15mG.mat;
load ./DANPOS_nucPositions_BY_60mG.mat;

gene_datasets = cell(5,1);
gene_datasets{1}=harry_annotations;
gene_datasets{2}=RPG_geneCoords;
gene_datasets{3}=OX_geneCoords;
gene_datasets{4}=RB_geneCoords;
gene_datasets{5}=RC_geneCoords;

dataset_labels = {'All genes','RPGs','OX genes','RB genes','RC genes'};

clear harry_annotations;
clear RPG_geneCoords;
clear OX_geneCoords;
clear RB_geneCoords;
clear RC_geneCoords;

%%
TSS_start=-20;

nuc_positions=cell(length(gene_datasets),1);
relative_shifts=cell(length(gene_datasets),1);
num_genes=zeros(length(gene_datasets),1);
for d=1:length(num_genes)
    for chr=1:size(gene_datasets{d},1)
        for strand=1:2
            num_genes(d)=num_genes(d)+size(gene_datasets{d}{chr,strand},1);
        end
    end
end

for d=1:length(gene_datasets)
    nuc_positions{d} = nan(3,4,num_genes(d));
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
                nuc_positions{d}(1,3,k)=min(idx_D((idx_D>=TSS_start+150)&(idx_D<=TSS_start+250)));
            end
            try
                nuc_positions{d}(1,4,k)=min(idx_D((idx_D>=TSS_start+300)&(idx_D<=TSS_start+400)));
            end
            
            try
                nuc_positions{d}(2,2,k)=min(idx_15mG((idx_15mG>=TSS_start)&(idx_15mG<=TSS_start+100)));
            end
            try
                nuc_positions{d}(2,3,k)=min(idx_15mG((idx_15mG>=TSS_start+150)&(idx_15mG<=TSS_start+250)));
            end
            try
                nuc_positions{d}(2,4,k)=min(idx_15mG((idx_15mG>=TSS_start+300)&(idx_15mG<=TSS_start+400)));
            end
            
            try
                nuc_positions{d}(3,2,k)=min(idx_60mG((idx_60mG>=TSS_start)&(idx_60mG<=TSS_start+100)));
            end
            try
                nuc_positions{d}(3,3,k)=min(idx_60mG((idx_60mG>=TSS_start+150)&(idx_60mG<=TSS_start+250)));
            end
            try
                nuc_positions{d}(3,4,k)=min(idx_60mG((idx_60mG>=TSS_start+300)&(idx_60mG<=TSS_start+400)));
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
                nuc_positions{d}(1,3,k)=min(idx_D((idx_D>=TSS_start+150)&(idx_D<=TSS_start+250)));
            end
            try
                nuc_positions{d}(1,4,k)=min(idx_D((idx_D>=TSS_start+300)&(idx_D<=TSS_start+400)));
            end
            
            try
                nuc_positions{d}(2,2,k)=min(idx_15mG((idx_15mG>=TSS_start)&(idx_15mG<=TSS_start+100)));
            end
            try
                nuc_positions{d}(2,3,k)=min(idx_15mG((idx_15mG>=TSS_start+150)&(idx_15mG<=TSS_start+250)));
            end
            try
                nuc_positions{d}(2,4,k)=min(idx_15mG((idx_15mG>=TSS_start+300)&(idx_15mG<=TSS_start+400)));
            end
            
            try
                nuc_positions{d}(3,2,k)=min(idx_60mG((idx_60mG>=TSS_start)&(idx_60mG<=TSS_start+100)));
            end
            try
                nuc_positions{d}(3,3,k)=min(idx_60mG((idx_60mG>=TSS_start+150)&(idx_60mG<=TSS_start+250)));
            end
            try
                nuc_positions{d}(3,4,k)=min(idx_60mG((idx_60mG>=TSS_start+300)&(idx_60mG<=TSS_start+400)));
            end
            
            k=k+1;
        end
    end
    
    relative_shifts{d}=nan(2,4,num_genes(d));
    p_vals = nan(2,4);
    for nuc=1:4
        for t=1:2
            relative_shifts{d}(t,nuc,:)=-nuc_positions{d}(1,nuc,:)+nuc_positions{d}(t+1,nuc,:);
            [~,p_vals(t,nuc)] = ttest(relative_shifts{d}(t,nuc,:));
        end
    end
    
    
    
    fig=figure('position',[0 0 1080 800]);
    ax1 = subplot('position',[0.01 0.05 0.23 0.8]);
    boxplot([reshape(relative_shifts{d}(1,1,:),1,size(relative_shifts{d},3));reshape(relative_shifts{d}(2,1,:),1,size(relative_shifts{d},3))]','orientation','horizontal');
    hold on;text(-20,1.3,['p = ',num2str(p_vals(1,1))],'fontsize',16);text(-20,2.3,['p = ',num2str(p_vals(2,1))],'fontsize',16);
    xlim([-50 50]);
    title('Minus One Nucleosome','fontsize',16);
    ax=gca;
    ax.YTickLabels = {'15mG','60mG'};
    ax.FontSize = 16;
    ax.YTickLabelRotation=90;
    
    ax2 = subplot('position',[0.26 0.05 0.23 0.8]);
    boxplot([reshape(relative_shifts{d}(1,2,:),1,size(relative_shifts{d},3));reshape(relative_shifts{d}(2,2,:),1,size(relative_shifts{d},3))]','orientation','horizontal');
    hold on;text(-20,1.3,['p = ',num2str(p_vals(1,2))],'fontsize',16);text(-20,2.3,['p = ',num2str(p_vals(2,2))],'fontsize',16);
    xlim([-50 50]);
    title('Plus One Nucleosome','fontsize',16);
    ax=gca;
    ax.YTickLabels = {'15mG','60mG'};
    ax.FontSize = 16;
    ax.YTickLabelRotation=90;
    
    subplot('position',[0.51 0.05 0.23 0.8]);
    boxplot([reshape(relative_shifts{d}(1,3,:),1,size(relative_shifts{d},3));reshape(relative_shifts{d}(2,3,:),1,size(relative_shifts{d},3))]','orientation','horizontal');
    hold on;text(-20,1.3,['p = ',num2str(p_vals(1,3))],'fontsize',16);text(-20,2.3,['p = ',num2str(p_vals(2,3))],'fontsize',16);
    xlim([-50 50]);
    title('Plus Two Nucleosome','fontsize',16);
    ax=gca;
    ax.YTickLabels = {'15mG','60mG'};
    ax.FontSize = 16;
    ax.YTickLabelRotation=90;
    
    subplot('position',[0.76 0.05 0.23 0.8]);
    boxplot([reshape(relative_shifts{d}(1,4,:),1,size(relative_shifts{d},3));reshape(relative_shifts{d}(2,4,:),1,size(relative_shifts{d},3))]','orientation','horizontal');
    hold on;text(-20,1.3,['p = ',num2str(p_vals(1,4))],'fontsize',16);text(-20,2.3,['p = ',num2str(p_vals(2,4))],'fontsize',16);
    xlim([-50 50]);
    title('Plus Three Nucleosome','fontsize',16);
    ax=gca;
    ax.YTickLabels = {'15mG','60mG'};
    ax.FontSize = 16;
    ax.YTickLabelRotation=90;
    
    stitle = suptitle(strcat(dataset_labels{d},', n = ',num2str(num_genes(d))));
    stitle.FontSize=24;
    
%     saveas(fig,strcat('./figures_20180216/',dataset_labels{d},'_nucPositions.fig'));
end
