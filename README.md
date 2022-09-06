# HostMicrobiome
Disturbances in microbial skin recolonization and cutaneous immune response following allogeneic stem cell transfer

General information:
This tool was developed to determine the co-localization of two cell types in microscopic images. The method was adapted from the spatial arrangement analysis-function of the DAIME software (https://www.microbial-ecology.net/daime). However, this script uses the centroid position of the cell populations of interest as an input, while the DAIME software requires the original microscopic images. Therefore, this script can be used with any image software as long as X and Y coordinates can be exported, which is a very common feature found in e.g. ImageJ or the StrataQuest image cytometry program (V6.01.167; TissueGnostics GmbH). 

In the original approach, the creation of linear dipoles results in a probability P(r) of the two ends (poles) to hit or miss the two populations in question and the sum of hits (Hr) is divided through the total number of hits and misses (Tr). In addition, the population densities (D1 and D2) are also considered and are calculated by dividing the total area of one population (A1 or A2) through the total area of the analysis region (At). Therefore, the DAIME software calculates the spatial relationships between two populations by: 
g(r)=P(r)/(2*D1*D2)
Where: P(r) = H(r) / T(r), D1 = A1 / At and D2 = A2 / At
Due to this difference in input data, we adjusted the original formula by calculating the total number of hits at range r as the difference between cumulative hits at or below range r and at or below range r-1 (H(r) = H_cumm(r)-H_cumm(r-1)). The total number of hits or misses at range r is given by the product between the number of possible dipoles at range r (the circumference of a circle of radius r: 2* ?*r) and the total area (At) (T(r)= 2* ?*r*At). In addition, as all our input objects are just the size of one pixel (= centroid) the population densities (D1 and D2) were calculated using the number instead of the area of one population (N1 or N2) through the total image area (At): D1 = N1/At, D2 = N2/At) With these changes, our adapted formula is the following: 
g(r)=([H_cumm(r)-H_cumm(r-1)]*At])/([(2 * pi * r) * (2 * N1 * N2)])

This script was generated using Python 3.9.7.

Hardware requirements:
This script is rather computationally intensive; therefore, a good processor is suggested. We used an Intel(R) Xeon(R) CPU E5-1630 v4 @ 3.70GHz, with 16 GB of RAM.

Parameters:
The script has the following parameters:
- ROI: Region Of Interest in microscopic image for which co-localization is to be determined. If left empty, all ROIs will be processed by default. 
- thr: Threshold to score co-localization (default 0.05).
- n: Number of iterations to be done (default 1000)
- op: Path to output filename.
- ip: Path to input directory.
- np: Number of processes to run in parallel (default = Available CPU in running system)
- str_input_file: Input file
- str_1st_class: String name of first class (=1st cell type)
- str_2nd_class: String name of second class (= 2nd cell type)
- strUse_same_class_interaction: Same class (=cell type) interaction: "yes" or different classes: "no"
- avg_neigh_um: When generating Hit values at range R, consider the average hits for all [R-avg_neigh_um, R+avg_neigh_um] ranges (default = 1.0)

Input file:
The first line of the file should include the following headers (Note that there is a tab character at the beginning of the string):
- ROI: Sample ID
- event_label: ID of the event (=cell within population)
- X_px: X coordinates of the centroid position in pixels
- Y_px: Y coordinates of the centroid position in pixels
- X_mm: X coordinates of the centroid position in millimeters
- Y_mm: Y coordinates of the centroid position in millimeters
- AOI: Size of the Area Of Interest (e.g. in square millimeters)
- Population: cell type of interest
- Region: ID of ROI (e.g. different tissue compartments)
- ROIc: Sample ID with underline characters

Next lines describe the events - one event per line (Note that there is a row index at the beginning of the string line).

Test files:
We have also included two test files: (1) One file with a lot of co-localization between the cells and (2) a second file with little to no interaction between the cells.

If you have any questions regarding the script, please contact: robert.nica@tissuegnostics.com

