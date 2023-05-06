rm -rf 0.*
rm -rf [1-9]*
rm -f *Transport.out
rm -rf processor*/0.*
rm -rf processor*/[1-9]*

find processor*/0 ! -name 'U' ! -name 'p' ! -name 'phi' -type f -exec rm -f {} +

