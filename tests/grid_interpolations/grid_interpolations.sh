make clean
make
rm output/*.smet
rm output/grids/*
# ./grid_interpolations 2008-12-01T00:00 2009-01-31T23:00
./grid_interpolations 2008-12-01T00:00 2008-12-05T00:00
