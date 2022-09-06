import os, random, argparse, sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import shapely.geometry as geometry
import argparse
from pathlib import Path
import time
import datetime
from datetime import datetime
import multiprocessing as mp
import math


# argument parsing
parser = argparse.ArgumentParser(description = "Simulate cells for fixed bac coordintates.")
parser.add_argument("--ROI", help = "Region of interest on picture for which simulation should be done. If empty, process all ROIs", default = "")
parser.add_argument("--thr", help = "Threshold to score interaction. Default 0.05.", default = 0.15)
parser.add_argument("--n", help = "Number of iterations to be done. Default = 1000.", default = 1000)
parser.add_argument("--op", help = "Path to output filename. Default = f:/InterScore/images/sim_res", default = "f:/InterScore/images/sim_res")
parser.add_argument("--ip", help = "Path to input directory. Default = f:/InterScore/images", default = "f:/InterScore/images")
parser.add_argument("--np", help = "Number of process to run in Parallel. Default = Available CPU in running system", default = mp.cpu_count() )
parser.add_argument("--str_input_file", help = "input file", default = "" )
parser.add_argument("--str_1st_class", help = "first class", default = "" )
parser.add_argument("--str_2nd_class", help = "second class", default = "" )
parser.add_argument("--strUse_same_class_interaction", help = "same class interaction (yes) or different classes (no)", default = "" )
parser.add_argument("--avg_neigh_um", help = "When generating Hit values at range R, consider the average hits for all [R-avg_neigh_um, R+avg_neigh_um] ranges", default = 1.0 )
#parser.add_argument("--strMethod_Standard", help = "same class interaction (yes) or different classes (no)", default = "" )
#parser.add_argument("--strMethod_DAIME", help = "same class interaction (yes) or different classes (no)", default = "" )

args = parser.parse_args()

# Get args
ROI, thr, n, op, ip, cpu_count1, str_input_file, str_1st_class, str_2nd_class, strUse_same_class_interaction, avg_neigh_um = args.ROI, args.thr, int(args.n), args.op, args.ip, int(args.np), args.str_input_file, args.str_1st_class, args.str_2nd_class, args.strUse_same_class_interaction, args.avg_neigh_um

if not os.path.isdir(op):
    os.makedirs(op)

os.chdir(ip)


# set input params for input files if params not set up
if ((len(str_input_file) == 0) or (len(str_1st_class) == 0) or (len(str_2nd_class) == 0) or (len(strUse_same_class_interaction) == 0)):
    #str_input_file  = "images_080120.txt"
    #str_1st_class = "bacs" # script uses the first class positions 
    #str_2nd_class = "leukos" # script uses the second class number for simulation
    #strUse_same_class_interaction = "no" # "yes" or "no"

    #str_input_file  = "PosCo_CD45pos.txt"
    #str_1st_class = "leukos" # script uses the first class positions 
    #str_2nd_class = "leukos" # script uses the second class number for simulation
    #strUse_same_class_interaction = "yes" # "yes" or "no"

    #str_input_file  = "NegCo_CD45neg.txt"
    #str_1st_class = "CD45_negative" # script uses the first class positions 
    #str_2nd_class = "CD45_negative" # script uses the second class number for simulation
    #strUse_same_class_interaction = "yes" # "yes" or "no"

    #str_input_file  = "STAPHY.txt"
    #str_1st_class = "bacs" # script uses the first class positions 
    #str_2nd_class = "leukos" # script uses the second class number for simulation
    #strUse_same_class_interaction = "no" # "yes" or "no"

    #str_input_file  = "STAPHY SerialBiop for Simulation_2.txt"
    #str_1st_class = "bacs" # script uses the first class positions 
    #str_2nd_class = "leukos" # script uses the second class number for simulation
    #strUse_same_class_interaction = "no" # "yes" or "no"

    #str_input_file  = "STAPHY aGVHD for Simulation_2.txt"
    #str_1st_class = "bacs" # script uses the first class positions 
    #str_2nd_class = "leukos" # script uses the second class number for simulation
    #strUse_same_class_interaction = "no" # "yes" or "no"

    #str_input_file  = "EUB338 GvHD- noGvHD for Simulation.txt"
    #str_1st_class = "bacs" # script uses the first class positions 
    #str_2nd_class = "leukos" # script uses the second class number for simulation
    #strUse_same_class_interaction = "no" # "yes" or "no"

    str_input_file  = "images_080120_fixed_2.txt"
    str_1st_class = "bacs" # script uses the first class positions 
    str_2nd_class = "leukos" # script uses the second class number for simulation
    strUse_same_class_interaction = "no" # "yes" or "no"

## set default method to run (if script params did no specify it)
#if (len(strMethod_Standard) == 0):
#    strMethod_Standard = "no"
#if (len(strMethod_DAIME) == 0):
#    strMethod_DAIME = "yes"


strMethod_DAIME = "yes"

sheets = pd.read_csv(str_input_file , sep = "\t").iloc[:, 1:]

str_in_file_name = Path(str_input_file).stem
op = op + "/" + str_in_file_name

if not os.path.isdir(op):
    os.makedirs(op)

# arrROIs_unique = sheets.ROIc.unique()
# print (sheets.ROIc)

## below are function defs

#get polygon based on input points in the testsheet
def get_polygon(testsheet):

    # testsheet = sheets.loc[(sheets.ROIc == ROI) & (sheets.population != "apcs")]

    # #split up testsheet
    # testsheet_1st_class = testsheet.loc[testsheet.population == str_first_class]
    # testsheet_2nd_class = testsheet.loc[testsheet.population == str_second_class]

    # #get number of bacs
    # nbacs = testsheet_bacs.shape[0]
    # nleukos = testsheet_leukos.shape[0]

    if(testsheet.shape[0]<5):
        return None

    #draw polygon around points
    point_collection = geometry.MultiPoint(list(zip(testsheet.X_mm, testsheet.Y_mm)))

    hull = point_collection.convex_hull

    #for buffering
    #increase hull size if it is too small
    true_area = testsheet.AOI.unique()[0]

    buffer = hull.area/100
    while hull.buffer(buffer).area < true_area:
        buffer += hull.area/200

    #create buffered hull
    hull_buffered = hull.buffer(buffer)

    return(hull_buffered)

# get random point inside polygon
def get_random_point_in_polygon(poly):
     minx, miny, maxx, maxy = poly.bounds
     while True:
         p = geometry.Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
         if poly.contains(p):
             return p

# generate a number of random cells (the same number as the number of n_2nd_class(bacs))
def get_random_cells(poly, n_2nd_class):

    rcells = []
    for _ in range(n_2nd_class):
        cell = get_random_point_in_polygon(poly)
        rcells.append(cell.coords[0])

    return np.vstack(rcells)

# run daime method
# see attached "daime formula.txt" for formula
def DAIME_score(testsheet_1st_class, testsheet_2nd_class, strUse_same_class_interaction, fPxSize, fMaxDist_radius_um, fStep_radius_um, area_ROI_px, avg_neigh_um):

    # len of output arrays
    fStep_radius_mm = fStep_radius_um / 1000
    nLen_px = int(round(fMaxDist_radius_um / fPxSize) + 2);
    nLen_um = int(round(fMaxDist_radius_um / fStep_radius_um) + 1);

    avg_neigh_px = avg_neigh_um / fPxSize
    avg_neigh_px = math.ceil(avg_neigh_px)
    avg_neigh_px = max(0, avg_neigh_px)

    # allocate output arrays
    arrRadius_px    = np.zeros(nLen_px, float)
    arrRadius_px_sq = np.zeros(nLen_px, float)
    arrH_cumm_px    = np.zeros(nLen_px, float)
    arrH_px         = np.zeros(nLen_px, float)
    arrG_px         = np.zeros(nLen_px, float)
    arrHAvgRange_px = np.zeros(nLen_px, float)
    # arrG_cumm_px    = np.zeros(nLen_px, float)
    N1 = testsheet_1st_class.shape[0]
    N2 = testsheet_2nd_class.shape[0]

    # init radius (in px)
    for k in range(1, nLen_px):
        arrRadius_px[k] = k
        arrRadius_px_sq[k] = arrRadius_px[k] * arrRadius_px[k]

    #arrRadius_um    = np.zeros(nLen_um, float)
    #arrRadius_mm    = np.zeros(nLen_um, float)
    #arrH_cumm_um    = np.zeros(nLen_um, float)
    #for k in range(1, nLen_um):
    #    arrRadius_um[k] = arrRadius_um[k - 1] + fStep_radius_um
    #    arrRadius_mm[k] = arrRadius_mm[k - 1] + fStep_radius_mm
    #    #arrRadius_mm_sq[k] = arrRadius_mm[k] * arrRadius_mm[k]
    #dist_curr = 0
    #dist_curr_sq = 0;
    #for x1, y1 in zip(testsheet_1st_class.X_mm, testsheet_1st_class.Y_mm):
    #    for x2, y2 in zip(testsheet_2nd_class.X_mm, testsheet_2nd_class.Y_mm):
    #        if ( (x1 != x2) or (y1 != y2) ):
    #            dist_curr = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    #            # dist_curr_sq = (x2 - x1)**2 + (y2 - y1)**2 # raise to square optimization
    #            for k in range(nLen_um):
    #                if dist_curr < arrRadius_mm[k]:
    #                # if dist_curr_sq < arrRadius_px_sq[k]: # raise to square optimization
    #                    arrH_cumm_um[k] = arrH_cumm_um[k] + 1

    # compute cumulative no of hits(interactions) at each range
    # (each object is assumed to be a single point - area is 1 square px)
    dist_curr = 0
    dist_curr_sq = 0;
    for x1, y1 in zip(testsheet_1st_class.X_px, testsheet_1st_class.Y_px):
        for x2, y2 in zip(testsheet_2nd_class.X_px, testsheet_2nd_class.Y_px):
            if ( (x1 != x2) or (y1 != y2) ):
                dist_curr = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
                # dist_curr_sq = (x2 - x1)**2 + (y2 - y1)**2 # raise to square optimization
                for k in range(nLen_px):
                    if dist_curr < arrRadius_px[k]:
                    # if dist_curr_sq < arrRadius_px_sq[k]: # raise to square optimization
                        arrH_cumm_px[k] = arrH_cumm_px[k] + 1

    # based on cumulative no of hits(interactions) compute the no of hits at each range
    for k in range(1, nLen_px):
        arrH_px[k] = arrH_cumm_px[k] - arrH_cumm_px[k - 1]

   # average hits with neigh range
    for k in range(0, nLen_px):
        arrHAvgRange_px[k] = 0
        nRangesCount = 0
        for j in range(-avg_neigh_px, avg_neigh_px + 1):
            k_neigh = k + j
            if ((k_neigh >= 0) and (k_neigh < nLen_px) ):
                arrHAvgRange_px[k] = arrHAvgRange_px[k] + arrH_px[k_neigh]
                nRangesCount = nRangesCount + 1
        if (nRangesCount > 0):
            arrHAvgRange_px[k] = arrHAvgRange_px[k] / nRangesCount

    # generate DAIME formula score in arrG_px
    val_mul_div = 2 * math.pi * (2 * N1 * N2)
    if val_mul_div != 0:
        for k in range(1, nLen_px):
            if arrRadius_px[k] > 0:
                # arrG[k] = (arrH[k] * area_ROI) / ((2 * pi * arrRadius[k]) * (2 * N1 * N2))
                # arrG_px[k]      = (arrH_px[k]              * area_ROI_px) / (val_mul_div * arrRadius_px[k])
                arrG_px[k]      = (arrHAvgRange_px[k]      * area_ROI_px) / (val_mul_div * arrRadius_px[k])
                # arrG_cumm_px[k] = (arrH_cumm_px[k] * area_ROI_px) / (val_mul_div * arrRadius_px[k])

    # allocate micrometer (output) arrays
    arrRadius_um = np.zeros(nLen_um, float)
    arrH_cumm_um = np.zeros(nLen_um, float)
    arrH_um      = np.zeros(nLen_um, float)
    arrG_um      = np.zeros(nLen_um, float)
    # arrG_cumm_um = np.zeros(nLen_um, float)

    # init radius as micrometer
    for k in range(1, nLen_um):
        arrRadius_um[k] = arrRadius_um[k - 1] + fStep_radius_um

    # convert from pixels to micrometers
    for k in range(1, nLen_um):
        val_um = arrRadius_um[k]
        val_px = val_um / fPxSize
        val_px = round(val_px)
        val_px = max(0, val_px)
        val_px = int(min(val_px, nLen_px - 1))
        arrH_cumm_um[k] = arrH_cumm_px[val_px]
        # arrH_um[k]      = arrH_px[val_px]
        arrH_um[k]      = arrHAvgRange_px[val_px]
        # arrG_cumm_um[k] = arrG_cumm_px[val_px]
        arrG_um[k]      = arrG_px[val_px]

    return arrRadius_um, arrH_cumm_um, arrH_um, arrG_um # , arrG_cumm_um

#get distances of bacs to cells (or type 1 to type 2)
#note X and Y coordinates of cells are passed to the function,
#testsheet bacs is defined above since bac coordinates are fixed
def dist_c2b(X_cells, Y_cells, thr, testsheet_1st_class, strUse_same_class_interaction):
    """
    Function that calculates the distance for each cell to bacteria (or generic type 1 to type 2).

    Parameters:
    -X_cells, Y_cells: X and Y coordinates for cells
    -thr: Threshold to identify interactors
    -testsheet_1st_class: Table for with coordinates.

    Returns:
    4-element tuple.
    1. fraction of interacting cells
    2. mean minimum distance
    3. std. min. distance
    4. median min. distance
    5. count of total interactions
    """

    #get number of bacs
    n_1st_class = testsheet_1st_class.shape[0]

    #list that holds whether or not a point falls with threshold
    ias = []

    #list that holds the minimum distances
    min_dists = []

    d_fraction_of_interacting_cells = 0
    d_mean_minimum_distance = 0
    d_std_min_distance = 0
    d_median_min_distance = 0
    d_total_interactions = 0

    if strUse_same_class_interaction == "yes":
        #calculate distance for all cells except self
        for x, y in zip(X_cells, Y_cells):
            dists = []
            ias_curr = False
            for x2, y2 in zip(testsheet_1st_class.X_mm, testsheet_1st_class.Y_mm):
                if ( (x != x2) or (y != y2) ):
                    dist_curr = np.sqrt((x - x2)**2 + (y - y2)**2)
                    dists.append(dist_curr)
                    if (dist_curr < thr):
                        ias_curr = True
                        d_total_interactions = d_total_interactions + 1
            ias.append(ias_curr)
            # ias.append(True) if (dists < thr).any() else ias.append(False)
            if (len(dists) > 0):
                min_dists.append(min(dists))

    else:
        n_count_interactions_curr = 0
        for x, y in zip(X_cells, Y_cells):
            #calculate distance for all cells
            dists = np.sqrt((x - testsheet_1st_class.X_mm)**2 + (y - testsheet_1st_class.Y_mm)**2)
            ias.append(True) if (dists < thr).any() else ias.append(False)
            dists_below_thr = (dists < thr);
            n_count_interactions_curr = np.sum(dists_below_thr)
            d_total_interactions = d_total_interactions + n_count_interactions_curr
            min_dists.append(min(dists))
            #mindist = np.min(dists)
            #print(mindist)
            #ias.append(mindist)

    if (len(ias) > 0):
        d_fraction_of_interacting_cells = np.sum(ias)/len(ias)
    if (len(min_dists) > 0):
        d_mean_minimum_distance = np.mean(min_dists)
        d_std_min_distance = np.std(min_dists)
        d_median_min_distance = np.median(min_dists)

    return d_fraction_of_interacting_cells, d_mean_minimum_distance, d_std_min_distance, d_median_min_distance, d_total_interactions

def main():

    ts = time.time()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print('\nStarted at: ', dt_string,'\n')

    # main part of script:

    # set the ROIs where measurements are generated (an unique ROI of all ROIs)
    arrROIs_unique = []
    if (len(ROI) > 0):
        arrROIs_unique = np.array([ROI])
    else:
        arrROIs_unique = sheets.ROIc.unique()

    # init dataframe
    output_df_allROIs = pd.DataFrame(columns = ['idx_ROI', 'ROI_Name', 'NoSimulations', 'avg_frac_obs', 'avg_mean_obs', 'avg_std_obs', 'avg_median_obs', 'avg_frac_sim', 'avg_mean_sim', 'avg_std_sim', 'avg_median_sim', 'pseudop_frac_ratio', 'pseudop_mean_ratio', 'pseudop_median_ratio', 'n_1st_class', 'n_2nd_class', 'area_roi', 'str_organ', 'avg_count_inter_sim', 'avg_count_inter_obs'])

    # get unique ROIs
    print ("Repeating apparitions of all ROIs = ", sheets.ROIc.size)
    print ("Unique ROIs = ", arrROIs_unique.size)
    # print (arrROIs_unique)
    nLen_unique = len(arrROIs_unique)

    # run DAIME-like adapted formula / method
    if strMethod_DAIME == "yes":
        # for each unique ROI in the input file

        output_df_DAIME_G = {}
        output_df_DAIME_G_cumm = {}
        output_df_DAIME_H = {}
        # output_df_DAIME_H_cumm = {}
        strRange_curr = []

        idx_roi = 0;
        idx_roi_end = arrROIs_unique.size - 1
        lstRanges = [];
        nGenerate_LstRanges = 0
        for idx_roi in range(0, nLen_unique):
            #get testsheet
            src_ROI=arrROIs_unique[idx_roi]
            testsheet = sheets.loc[(sheets.ROIc == src_ROI) & (sheets.population != "apcs")]
            #split up testsheet
            testsheet_1st_class = testsheet.loc[testsheet.population == str_1st_class]
            testsheet_2nd_class = testsheet.loc[testsheet.population == str_2nd_class]

            #init parameters from testsheet
            fPxSize = 1
            area_roi_mm = 0
            str_organ = "missing"
            if testsheet.shape[0] > 0:
                area_roi_mm = testsheet['AOI'].astype(float).mean()
                str_organ = testsheet.iloc[0]['region']
                # the fPxSize_x and fPxSize_y formula for each point should give the same value
                pt_x_mm = testsheet['X_mm'].astype(float).mean()
                pt_y_mm = testsheet['Y_mm'].astype(float).mean()
                pt_x_px = testsheet['X_px'].astype(float).mean()
                pt_y_px = testsheet['Y_px'].astype(float).mean()
                # compute the pixel size
                fPxSize_x = (pt_x_mm * 1000) / pt_x_px
                fPxSize_y = (pt_y_mm * 1000) / pt_y_px
                fPxSize = (fPxSize_x + fPxSize_y) / 2

            fMaxDist_radius_um = thr * 1000
            fStep_radius_um = 1
            area_roi_um = area_roi_mm * 1000 * 1000
            # A_um = A_px * px_sz^2
            # A_px = A_um / px_sz^2
            area_ROI_px = area_roi_um / (fPxSize * fPxSize)

            arrRadius, arrH_cumm, arrH, arrG = DAIME_score(testsheet_1st_class, testsheet_2nd_class, strUse_same_class_interaction, fPxSize, fMaxDist_radius_um, fStep_radius_um, area_ROI_px, avg_neigh_um)
            nLenArrays = len(arrRadius)

            # generate output dfs
            if nGenerate_LstRanges == 0:
                lstRanges.append("region")
                lstRanges.append("N1")
                lstRanges.append("N2")
                lstRanges.append("area_mm")
                lstRanges.append("thr_mm")
                for idx_range in range(0, nLenArrays):
                    strRange_curr = "range_" + str(arrRadius[idx_range]) + "um";
                    lstRanges.append(strRange_curr)
                output_df_DAIME_G      = pd.DataFrame(columns = lstRanges, index = range(nLen_unique))
                # output_df_DAIME_G_cumm = pd.DataFrame(columns = lstRanges, index = range(nLen_unique))
                output_df_DAIME_H      = pd.DataFrame(columns = lstRanges, index = range(nLen_unique))
                output_df_DAIME_H_cumm = pd.DataFrame(columns = lstRanges, index = range(nLen_unique))
                nGenerate_LstRanges = 1

            # init current src_ROI
            idx_col_aux = 0
            output_df_DAIME_G.iat     [idx_roi, idx_col_aux] = src_ROI
            # output_df_DAIME_G_cumm.iat[idx_roi, idx_col_aux] = src_ROI
            output_df_DAIME_H.iat     [idx_roi, idx_col_aux] = src_ROI
            output_df_DAIME_H_cumm.iat[idx_roi, idx_col_aux] = src_ROI

            # init current type 1 no
            idx_col_aux = 1
            output_df_DAIME_G.iat     [idx_roi, idx_col_aux] = str(testsheet_1st_class.shape[0])
            # output_df_DAIME_G_cumm.iat[idx_roi, idx_col_aux] = str(testsheet_1st_class.shape[0])
            output_df_DAIME_H.iat     [idx_roi, idx_col_aux] = str(testsheet_1st_class.shape[0])
            output_df_DAIME_H_cumm.iat[idx_roi, idx_col_aux] = str(testsheet_1st_class.shape[0])

            # init current type 2 no
            idx_col_aux = 2
            output_df_DAIME_G.iat     [idx_roi, idx_col_aux] = str(testsheet_2nd_class.shape[0])
            # output_df_DAIME_G_cumm.iat[idx_roi, idx_col_aux] = str(testsheet_2nd_class.shape[0])
            output_df_DAIME_H.iat     [idx_roi, idx_col_aux] = str(testsheet_2nd_class.shape[0])
            output_df_DAIME_H_cumm.iat[idx_roi, idx_col_aux] = str(testsheet_2nd_class.shape[0])

            # init current ROI area (in mm)
            idx_col_aux = 3
            output_df_DAIME_G.iat     [idx_roi, idx_col_aux] = str(area_roi_mm)
            # output_df_DAIME_G_cumm.iat[idx_roi, idx_col_aux] = str(area_roi_mm)
            output_df_DAIME_H.iat     [idx_roi, idx_col_aux] = str(area_roi_mm)
            output_df_DAIME_H_cumm.iat[idx_roi, idx_col_aux] = str(area_roi_mm)

            # init current threshold
            idx_col_aux = 4
            output_df_DAIME_G.iat     [idx_roi, idx_col_aux] = str(thr)
            # output_df_DAIME_G_cumm.iat[idx_roi, idx_col_aux] = str(thr)
            output_df_DAIME_H.iat     [idx_roi, idx_col_aux] = str(thr)
            output_df_DAIME_H_cumm.iat[idx_roi, idx_col_aux] = str(thr)

            # init current results for all ranges up to threshold
            for idx_range in range(0, nLenArrays):
                # output_df_DAIME_G.at[idx_range, idx_roi] = arrG[idx_range]
                output_df_DAIME_G.iat     [idx_roi, idx_range + idx_col_aux + 1] = str(arrG     [idx_range])
                # output_df_DAIME_G_cumm.iat[idx_roi, idx_range + idx_col_aux + 1] = str(arrG_cumm[idx_range])
                output_df_DAIME_H.iat     [idx_roi, idx_range + idx_col_aux + 1] = str(arrH     [idx_range])
                output_df_DAIME_H_cumm.iat[idx_roi, idx_range + idx_col_aux + 1] = str(arrH_cumm[idx_range])

            print('ROI: ', src_ROI, idx_roi, '/', idx_roi_end, 'processed', '\n')

        #write result files:
        op_file_daime_g = os.path.join(op, Path(str_input_file).stem + "_DAIME_G_all_ROIs.txt")
        output_df_DAIME_G.to_csv(op_file_daime_g, sep = "\t")

        #op_file_daime_gcumm = os.path.join(op, Path(str_input_file).stem + "_DAIME_G_cumm_all_ROIs.txt")
        #output_df_DAIME_G_cumm.to_csv(op_file_daime_gcumm, sep = "\t")

        op_file_daime_h = os.path.join(op, Path(str_input_file).stem + "_DAIME_H_all_ROIs.txt")
        output_df_DAIME_H.to_csv(op_file_daime_h, sep = "\t")

        op_file_daime_hcumm = os.path.join(op, Path(str_input_file).stem + "_DAIME_H_cumm_all_ROIs.txt")
        output_df_DAIME_H_cumm.to_csv(op_file_daime_hcumm, sep = "\t")

    time_in_sec = time.time() - ts
    time_in_min = time.strftime("%H:%M:%S", time.gmtime(time_in_sec))
    print('\nTotal run time:', time_in_min,'\n')

if __name__ == '__main__':
    main()
