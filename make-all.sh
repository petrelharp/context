
for i in $(echo */makefile)
do
    make -C $(echo $i | sed -e "s/makefile//")
done
