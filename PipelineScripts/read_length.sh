
echo "RUN FOR MULTIPLE FILES AS FOLLOWS: "
echo "find . -name \"\*.gz\" -exec read_length.sh {} \\;"

echo $1
zless $1 | head -100000 | awk '{if ((NR%4)==0) printf "%s\n",$0; else printf "%s\t",$0;}' | awk '{print length($3)}' | sort | uniq -c 
