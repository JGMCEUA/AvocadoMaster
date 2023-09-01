!/bin/sh
for x in *; do tar -czvf $x.tar.gz $x; done
