#exp=xaema
#exp=xaemb
#exp=xaemc

exp=xaemd

#file=pd_???.nc
file=.nc

#for suffix in a b c d  e f g h i j k l m n o p q r s t u v w x y z 
for suffix in  a b c d  e f g h i j k l m n o p q r s t u v
do
 
mkdir $exp$suffix                                                    
cd  $exp$suffix

echo "cd /data/brussel/vo/iceclim/vsc10103/ISM/aism/eot/R/output_emul_20geom/$exp$suffix" > msub.tmp
echo "get *$file" >> msub.tmp                                             
echo "exit" >> msub.tmp
    
sftp -b msub.tmp vsc
rm msub.tmp
cd ../


done
