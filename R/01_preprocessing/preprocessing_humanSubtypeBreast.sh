# make directory for the data first time 
# mkdir outputs/filtered
# cd outputs/filtered

mkdir humanSubtypeBreast
cd humanSubtypeBreast  

mkdir outs

wget "https://zenodo.org/records/4739739/files/spatial.tar.gz?download=1" -O spatial.tar.gz
wget "https://zenodo.org/records/4739739/files/filtered_count_matrices.tar.gz?download=1" -O filtered_count_matrices.tar.gz

gunzip filtered_count_matrices.tar.gz
tar -xvf filtered_count_matrices.tar

gunzip spatial.tar.gz
tar -xvf spatial.tar

# clean up
rm filtered_count_matrices.tar
rm spatial.tar

# copy over the data from one sample to outs folder
cp filtered_count_matrices/CID4290_filtered_count_matrix outs -r
cp spatial/CID4290_spatial outs -r

# clean up
cd outs
mv CID4290_filtered_feature_bc_matrix filtered_feature_bc_matrix
mv CID4290_spatial spatial

# download meta data to remove artefacts
wget "https://zenodo.org/records/4739739/files/metadata.tar.gz?download=1" -O metadata.tar.gz
gunzip metadata.tar.gz
tar -xvf metadata.tar
rm metadata.tar
