'''
Due to issues with the past analysis script its difficult to gop directorys when evaulutaing the md data.
The purpose of this script is to recursively apply it to all the complex directorys and get any data from the md dirs.
It needs to do the following:
1.) Find the md directorys in the omplex directorys
2.) run md analysis
    a.) RMSD
        i.) Graphs?
    b.) Lysine accessibility
        i.) it would be a huge bonus if this works but it may need a lot of time to develop
    c.) compile MMGBSA data
        i.) comparitive data between no linker control would be nice
3.) Output data to output.txt, previously established from analysis.py
4.) Create a folder to bundle everything and zip it for use off the clusters
    a.) output.txt
    b.) graphs
    c.) mmgbsa.out
    d.) top5complexes.csv
'''