# GMX_Hbonds_plot

In the repository you will find two scripts in the scripts folder.
- hbonds.sh
- hbonds.py

**<a href="https://github.com/mangeshdamre/GMX_Hbonds_plot/blob/main/scripts/hbonds.sh" target="_blank">hbonds.sh</a>** will generate the data for the given selection in the script. One can update the selection for hbonds calculation in the script at the section "gmx make_ndx"
```sh
gmx make_ndx -f ../$GRO -o hbonds.ndx << EOL
a 2708-2792
name 10 A

a 5852-5936
name 11 B

a 8996-9080
name 12 C

a 12140-12224
name 13 D

a 15284-15368
name 14 E

a 18428-18512
name 15 F

q
EOL
```

The script is for protein assuming 6 protomers (chains=A..F). One can update the chains as per the protein of interest.
```sh
for chain in A B C D E F  
do
mkdir $chain
if [ ! -f $chain/hbnum$chain.xvg ]; then
	gmx hbond -f ../$XTC \
	-s ../../$TPR \
	-n hbonds.ndx \
	-num $chain/hbnum$chain.xvg \
	-g $chain/hbond$chain.log \
	-ac $chain/hbac$chain.xvg \
	-dist $chain/hbdist$chain.xvg \
	-hx $chain/hbhelix$chain.xvg \
	-hbn $chain/hbond$chain.ndx \
	-hbm $chain/hbmap$chain.xpm \
	-a 20 -r 0.3 <<EOL
$chain
$chain
EOL
fi
done
```
