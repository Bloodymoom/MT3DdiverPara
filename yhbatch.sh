#!/bin/bash

NODES=$1
PRO=$2
MAXIT=$3
TOLE=$4
TYPE=$5

cat > yhrun.sh << EOF
#!/bin/bash

yhrun -pthcp1 -N $NODES -n $PRO ./main -ksp_max_it $MAXIT -ksp_rtol $TOLE
EOF

yhbatch -pthcp1 -N $NODES -n $PRO -o log/log-$TYPE-$PRO-n$MAXIT.out yhrun.sh

