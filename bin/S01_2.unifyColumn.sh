#usage: sh $0 $1 $2 > outfile
#$1 and $2 are two file needed to output the common column
head -1 $1|xargs -n1 > unify.tmp.1
head -1 $2|xargs -n1 > unify.tmp.2

comm $unify.tmp.f1 $unify.tmp.f2 -1 -2


