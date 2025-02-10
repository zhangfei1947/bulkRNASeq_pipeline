import pandas as pd
import os,sys,csv


# 定义输入和输出目录
anno_file = sys.argv[1]
input_dir = sys.argv[2]
postfix = sys.argv[3]
output_dir = sys.argv[4]

# 读取基因注释文件
annotation_file = anno_file
annotation_df = pd.read_csv(annotation_file, sep='\t', quotechar='"')

# 遍历输入目录中的所有差异分析结果文件
for filename in os.listdir(input_dir):
    if filename.endswith(postfix):
        # 构建完整的文件路径
        diff_expr_file = os.path.join(input_dir, filename)
        
        # 读取差异分析结果文件，提取 geneID 和 padj 列
        diff_expr_df = pd.read_csv(diff_expr_file, sep='\t', usecols=['padj'])
        diff_expr_df = diff_expr_df.reset_index()
        diff_expr_df = diff_expr_df.rename(columns={'index': 'geneID'})
        
        # 合并两个数据框，基于 geneID
        merged_df = pd.merge(diff_expr_df, annotation_df, left_on='geneID', right_on='ID')
        
        # 选择所需的列，并重新排列顺序
        result_df = merged_df[['geneID', 'padj', 'Name', 'fullname', 'Ontology_term', 'Alias', 'Dbxref', 'gbunit']]
        
        # 构建输出文件路径
        output_file = os.path.join(output_dir, filename)
        
        # 将结果写入新的文件
        result_df.to_csv(output_file, sep='\t', quotechar='"', quoting=csv.QUOTE_ALL, index=False)

print(f"Annotated results have been saved to {output_dir}")
