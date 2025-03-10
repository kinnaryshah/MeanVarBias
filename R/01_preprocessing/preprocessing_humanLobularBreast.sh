# make directory for the data first time 
# mkdir outputs/raw
cd outputs/raw

mkdir humanLobularBreast
cd humanLobularBreast  

mkdir outs
cd outs

wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_raw_feature_bc_matrix.tar.gz
wget https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_spatial.tar.gz

gunzip Parent_Visium_Human_BreastCancer_raw_feature_bc_matrix.tar.gz
tar -xvf Parent_Visium_Human_BreastCancer_raw_feature_bc_matrix.tar

gunzip Parent_Visium_Human_BreastCancer_spatial.tar.gz
tar -xvf Parent_Visium_Human_BreastCancer_spatial.tar

# clean up
rm Parent_Visium_Human_BreastCancer_raw_feature_bc_matrix.tar
rm Parent_Visium_Human_BreastCancer_spatial.tar
