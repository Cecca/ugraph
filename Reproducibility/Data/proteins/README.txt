= Proteins datasets =

This folder contains datasets of protein-protein interactino networks.
Files `collins2007.txt`, `gavin2006_socioaffinities_rescaled.txt`,
`krogan2006_core.txt`, and `krogan2006_extended.txt` have been
downloaded from
http://www.paccanarolab.org/static_content/clusterone/cl1_datasets.zip

The TAP dataset has been downloaded from
http://tap.med.utoronto.ca/exttap/downloads.php and is the same as
`krogan2006_core.txt`, with a different node labelling.
The dataset is comprised of the following files
 - Krogan TAP uncertain graph:
   http://tap.med.utoronto.ca/exttap/downloads/TAP_core.txt, see file
   `TAP_core.txt`.
 - Krogan MIPS ground truth: derived with the
   `scripts/krogan_ground_truth.py` script from the csv version of the
   excel file that can be downloaded from
   http://tap.med.utoronto.ca/exttap/downloads/MIPS_annotations_for_Krogan-etal_Complexes.xls.
   See file `TAP_ground.txt`.
 - Krogan clustering using MCL: derived from
   http://tap.med.utoronto.ca/exttap/downloads/MCL_clusters.txt
   using the following sequence of python code
   
     fp = open("MCL_clusters.txt")
     orig = fp.readlines()[1:]
     fp.close()
     new_lines = [" ".join(l.split()[2:]) + "\n" for l in orig]
     ofp = open("Krogan-clusters.txt", "w")
     ofp.writelines(new_lines)
     ofp.close()

   See file `TAP_clusters.txt`.
   
 

