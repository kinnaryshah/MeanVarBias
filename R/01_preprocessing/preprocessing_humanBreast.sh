# make directory for the data first time 
# mkdir outputs/raw
# cd outputs/raw

mkdir humanBreast
cd humanBreast  

mkdir outs
cd outs

wget https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Human_Breast_Cancer/Visium_Human_Breast_Cancer_spatial.tar.gz

gunzip Visium_Human_Breast_Cancer_raw_feature_bc_matrix.tar.gz
tar -xvf Visium_Human_Breast_Cancer_raw_feature_bc_matrix.tar

gunzip Visium_Human_Breast_Cancer_spatial.tar.gz
tar -xvf Visium_Human_Breast_Cancer_spatial.tar
