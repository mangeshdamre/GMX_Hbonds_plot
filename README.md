# GMX_Hbonds_plot

In the repository you will find two scripts in the scripts folder.
- hbonds.sh
- hbonds.py

**<a href="https://github.com/mangeshdamre/GMX_Hbonds_plot/blob/main/scripts/hbonds.sh" target="_blank">hbonds.sh</a>**) will generate the data for the given selection in the script. One can update the selection for hbonds calculation in the script at the section "gmx make_ndx"
```
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
