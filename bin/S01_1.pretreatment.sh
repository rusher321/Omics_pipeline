### configuration
 input_csv=$1
 target_ID=$2
current_ID=$3
    prefix=$4
#ID_pfx means the STRING PART of the sample ID, which used to split the ID and sort according to the numeric part
	ID_pfx=$5

#lib
  idmap_pl=$6
replace_pl=$7
   trim_pl=$8
   sort_pl=$9

#Annotate the ID* steps if you have already unified the ID between intensity data and the phenotype sheet

# ID1 # build ID mapping (And it may need to be modified so as to fit your data)
perl $idmap_pl $target_ID $current_ID > $prefix.map

# 01 # remove quote and tranform csv to xls
sed 's/,/\t/g;s/"//g' $input_csv > $prefix.xls

# ID2 # transform ID to Blood Barcode
perl $replace_pl $prefix.map $prefix.xls h > $prefix.renumber.xls

# 02 # remove the redundancy row & label line
perl $trim_pl $prefix.map $prefix.renumber.xls in $prefix.renumber.trimed.xls h

sed -i '2d' $prefix.renumber.trimed.xls

# 03 # sort column ID
perl $sort_pl $prefix.renumber.trimed.xls $ID_pfx	# > $prefix.renumber.trimed.sorted.xls

exit
