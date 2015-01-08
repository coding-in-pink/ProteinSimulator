MI 214
Programming Project 3

Instructions: 
Run the md program with 'python md_em073094.py --iF input_file --kB --kN --nbCutoff --m --dt --n --out' (where all but the first argument are optional)

GRAPHS/EUC:
For the two .euc files, I created a separate method in my program that wrote the euclidean distances between the two points (either site1 or site2) for a specified step to a separate euc file. 
Because this file only corresponds to a single .rvc file, I ran this for each of the three .rvc files and then transferred the data to Matlab. 
In Matlab, I concatenated the three data vectors together, and plotted them to make the two graphs. I then transferred the concatenated data vectors back to the site1/site2.euc files so the data for the 3 .rvc files appeared side by side.
