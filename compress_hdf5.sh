
path=$1;
cd ${path}
for f in $(find -name "*.hdf5"); do
    echo ${f};
    h5repack -v -f GZIP=9  ${f} "compressed"${f}
done
