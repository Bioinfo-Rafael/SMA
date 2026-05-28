#First File
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287569
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173524
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167332
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219201
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242942
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208629
¥

obs: cell_id / gse_id / sample_id / condition / cell_type(無ければ 'unknown')を統一カラムに正規化、元のobsカラムも保持
var_names: 遺伝子シンボル優先（無ければENSEMBL ID） + gene_symbol_upper で大文字小文字差を吸収する比較用カラム
X: 生UMIカウントのCSR sparse
セルバーコードは {GSM}_{barcode} でグローバルにユニーク化