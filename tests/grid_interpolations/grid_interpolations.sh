make clean > /dev/null 2>&1
make > /dev/null
rm output/*.smet > /dev/null 2>&1
rm output/grids/* > /dev/null 2>&1
./grid_interpolations 2008-12-01T08:00 2008-12-05T00:00
