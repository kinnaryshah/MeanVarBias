cd outputs/raw

mkdir humanOvarian
cd humanOvarian

mkdir outs
cd outs

wget -O spatial.zip "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6506nnn/GSM6506111/suppl/GSM6506111_SP2_spatial.zip"
unzip spatial.zip
rm spatial.zip

mkdir filtered_feature_bc_matrix
cd filtered_feature_bc_matrix

wget -O barcodes.tsv.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6506nnn/GSM6506111/suppl/GSM6506111_SP2_barcodes.tsv.gz"
wget -O features.tsv.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6506nnn/GSM6506111/suppl/GSM6506111_SP2_features.tsv.gz"
wget -O matrix.mtx.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6506nnn/GSM6506111/suppl/GSM6506111_SP2_matrix.mtx.gz"


