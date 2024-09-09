cd outputs/raw

mkdir humanHPC_V12D07-335_D1
cd humanHPC_V12D07-335_D1

mkdir outs
cd outs

mkdir raw_feature_bc_matrix
cd raw_feature_bc_matrix

wget -O barcodes.tsv.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8226nnn/GSM8226202/suppl/GSM8226202_V10B01-085_A1_barcodes.tsv.gz"
wget -O features.tsv.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8226nnn/GSM8226202/suppl/GSM8226202_V10B01-085_A1_features.tsv.gz"
wget -O matrix.mtx.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8226nnn/GSM8226202/suppl/GSM8226202_V10B01-085_A1_matrix.mtx.gz"

cd ..
mkdir spatial
cd spatial

wget -O tissue_lowres_image.png.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8226nnn/GSM8226202/suppl/GSM8226202_V10B01-085_A1-tissue_lowres_image.png.gz"
gunzip tissue_lowres_image.png.gz

wget -O tissue_positions_list.csv.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8226nnn/GSM8226202/suppl/GSM8226202_V10B01-085_A1-tissue_positions_list.csv.gz"
gunzip tissue_positions_list.csv.gz

wget -O scalefactors_json.json.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8226nnn/GSM8226202/suppl/GSM8226202_V10B01-085_A1-scalefactors_json.json.gz"
gunzip scalefactors_json.json.gz

wget -O detected_tissue_image.jpg.gz "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8226nnn/GSM8226202/suppl/GSM8226202_V10B01-085_A1-detected_tissue_image.jpg.gz"
gunzip detected_tissue_image.jpg.gz
