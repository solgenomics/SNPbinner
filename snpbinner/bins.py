"""This script takes an output file from the crosspoints script (X.crosspoints.csv) and a minimum bin size and identifies a list of bins and their genotypes for each RIL. In order to do so, the crosspoints from all RILs are projected onto one representative chromosome. These crosspoints are then partitioned into groups where the distance between each consectutive crosspoint is less than the minimum bin size. Groups containing less than three crosspoint are combined into one representitive crosspoint located at the centroid of the group. If the group contains three or more crosspoints, the maximum number of breakpoints that fit inside the region bounded by the first and last crosspoint in each list is first determined by dividing the region size by the minimum bin size (rounded up). Then, clustering isperformed for all values of K from the maximum number of breakpoints to 1 using a modified 1D K-means algorithm which adjusts the centroid values such that they are greater than a minimum distance from eachother during the update step. Each of the produced clusterings are then scored using the average deviation for each cluster to the adjusted centroid. The adjusted centroids from the clustering with the lowest average deviation are then considered to be the representitive crosspoints for the group of crosspoints. The representative crosspoints determined for each group are then used as the bounds of the bins created. To genotype the bins for each RIL, the genotype which covers the most area inside of the bin fromthe original crosspoint data is used. The bin locations, bounds, and genotypes for each RIL are then output in CSV format."""

from collections import OrderedDict
import math

def parser(parser_add_func,name):
    '''Sets up an arguement parser for this module. Note the arguement names match those in the `run` function.'''
    p = parser_add_func(name,description=__doc__)
    p.add_argument("-i","--input",        metavar="PATH",  dest='input_path',   required=True,           help="Path to a crosspoints CSV.")
    p.add_argument("-o","--output",       metavar="PATH",  dest='output_path',  required=True,           help="Path for the output CSV.")
    p.add_argument("-l","--min-bin-size", metavar="INT",   dest='min_bin_size', required=True, type=int, help="Minimum size of a bin in basepairs. This defines the resolution of the binmap.")
    p.add_argument("-n","--binmap-id",    metavar="ID",    dest="binmap_id",    default=False, type=str, help="If a binmap ID is provided, a header row will be added and each column labeled with the given string.")
    return p

def run(input_path,output_path,min_bin_size,binmap_id):
    '''Runs the module'''

    #load all RIL crosspoint data and find the chromosome length
    line_cps = OrderedDict()
    chrom_len = 0
    with open(input_path) as input_file:
        for line in input_file:
            name,cps = line.split(",",1)
            cps = [item.strip() for item in cps.split(",") if item.strip()!=""]
            for i in range(0,len(cps),2):
                cps[i] = int(cps[i])
            line_cps[name] = cps
            if cps[-1] > chrom_len: chrom_len = cps[-1]

    #clears the output file and also ensure that it exists before continuing.
    with open(output_path,"w") as clear_output:
        pass

    # find all crosspoints and count how many times crosspoints occur at each location.
    cp_loc_count = {chrom_len:float('inf'),0:float('inf')}
    for cps in line_cps.values():
        for i in range(0,len(cps),2):
            if not cps[i] in cp_loc_count: cp_loc_count[cps[i]] = 0
            cp_loc_count[cps[i]] += 1
    all_cp = sorted(cp_loc_count.keys())
    print all_cp

    # partition all cp into groups where the distance between each consecutive cp is less than min_bin_size
    crosspoint_groups = [[all_cp[0]]]
    i = 1
    while i<len(all_cp):
        if all_cp[i]-crosspoint_groups[-1][-1] < min_bin_size:
            crosspoint_groups[-1].append(all_cp[i])
        else:
            crosspoint_groups.append([all_cp[i]])
        i+=1
    print crosspoint_groups

    bin_bounds = [] #Create a list for storing the bin bounds
    crosspoint_groups.reverse() #Moving backwards through the discovered groups

    #for each crosspoint group, predict bin bounderies using k-means
    for i,group in enumerate(crosspoint_groups):

        group_len = group[-1]-group[0]
        expanded_group = [cp for cp in group for num in range(cp_loc_count[cp] if cp_loc_count[cp]!=float('inf') else 1)]

        print "\n\n"+("="*80)
        print ("%s crosspoints, start=%s, end=%s"%(len(expanded_group),group[0],group[-1]))
        print "."*80

        #if there is only one cp in the group, simply use that as the one representitive boundary
        if group_len==0:
            print "\nN=1"
            bin_bounds.append(group[0])
            print "F", bin_bounds[-1]
            continue

        #Calculate the maximum number of new boundaries that could fit from the first to last crosspoint of the group
        max_new_cp = group_len//min_bin_size +1

        #if it is less than two, simply average the crosspoints for the best fitting boundary
        if max_new_cp<2:
            print "\nN=1"
            print "g",bin_bound_visualize(expanded_group,group[0],group[-1])
            bin_bounds.append(crosspoint_avg(cp_loc_count,group,chrom_len))
            print "f",bin_bound_visualize(bin_bounds[-1:],group[0],group[-1],aura=min_bin_size)
            continue

        #For each possible number of boundaries, predict the locations using k-means and the sum of the variance between each crosspoint to the closest boundary
        solution_list = []
        for cp_count in range(max_new_cp,0,-1):

            print "\n"+("N=%s"%cp_count)+"\n"
            print "g",bin_bound_visualize(expanded_group,group[0],group[-1])

            #initilize the k-means psuedo-centroids (they are not true centroids once adjusted for minimum distance) to be evenly spaced within the group
            start_cp_dist = group_len/float(cp_count)
            if min_bin_size>start_cp_dist:
                start_cp_dist = min_bin_size
            start_offset = group_len/2.0-(cp_count*start_cp_dist)/2.0+start_cp_dist/2
            km_points = [int((j*start_cp_dist)+group[0]+start_offset) for j in range(cp_count)]
            km_groups = [[] for k in range(len(km_points))]

            #assign each cp to the closest centroid, we are working in one dimension, so this can be done linearly (hence the upcoming while loop)
            nearest = 0
            for cp in group:
                while nearest+1<len(km_points) and abs(cp-km_points[nearest+1]) < abs(cp-km_points[nearest]):
                    nearest+=1
                km_groups[nearest].append(cp)
            print "u",bin_bound_visualize(km_points,group[0],group[-1],aura=min_bin_size)

            #Now that the cps have been assigned to km_groups (groups with a common closest centroid), perform k-means optimization!
            memo = set() #stores each visited state so that minima cycles can be detected
            adjustment_needed = True #remains true while we have not reached a final k-means result
            while adjustment_needed:
                change="non" #This string is for print output only.

                #recalculate the centroids
                km_points = [crosspoint_avg(cp_loc_count,k,chrom_len) for k in km_groups]

                #Checks to make sure there were no empty groups, if there were, it places the centroid of the empty group at the midpoint between the surrounding centroids.
                for k in range(1,len(km_points)-1):
                    if math.isnan(km_points[k]):
                        next_not_nan = k+1
                        while math.isnan(km_points[next_not_nan]):
                            next_not_nan+=1
                        km_points[k] = (km_points[k-1]+km_points[next_not_nan])//2

                #Simple(ish) lamda function to check for overlaps between minimum-sized bins around each centroid
                get_overlaps = lambda:[(overlap,k) for overlap,k in ((min_bin_size-abs(km_points[k]-km_points[k-1]),k) for k in range(1,len(km_points))) if overlap>0]
                overlap_adjustment_performed = False

                #create an id_tuple for current state of the 
                id_tuple = (tuple(len(g) for g in km_groups),tuple(km_points))
                if id_tuple in memo:
                    overlaps = get_overlaps()
                    while len(overlaps)>0:
                        for overlap,k in overlaps:
                            if overlap>0:
                                overlap_adjustment_performed=True
                                left_shift = overlap-(overlap//2)
                                right_shift = left_shift
                                if k==1 and left_shift>km_points[k-1]-group[0]:
                                    left_shift = km_points[k-1]-group[0]
                                    right_shift = overlap-left_shift
                                elif k==len(km_points)-1 and right_shift>group[-1]-km_points[k]:
                                    right_shift = group[-1]-km_points[k]
                                    left_shift = overlap-right_shift
                                km_points[k] += right_shift
                                km_points[k-1] -= left_shift
                        overlaps = get_overlaps()
                    #changed the centroids, so the ID changes
                    id_tuple = (tuple(len(g) for g in km_groups),tuple(km_points))

                #reasign each cp to the km_group belonging to each centroid! Again, is done in linear time.
                for k in range(len(km_groups)):
                    if k>0:   
                        while len(km_groups[k])>0 and abs(km_groups[k][0]-km_points[k-1]) < abs(km_groups[k][0]-km_points[k]):
                            km_groups[k-1].append(km_groups[k].pop(0))
                    if k<len(km_groups)-1:  
                        while len(km_groups[k])>0 and abs(km_groups[k][-1]-km_points[k+1]) < abs(km_groups[k][-1]-km_points[k]):
                            km_groups[k+1].insert(0,km_groups[k].pop())

                #if we have been here before, we are done!
                if id_tuple in memo:
                        adjustment_needed = False
                        print "Done."
                else:
                    memo.add(id_tuple)
                    if overlap_adjustment_performed: 
                        change="ovr"
                    else:
                        change="adj"
                print change,bin_bound_visualize(km_points,group[0],group[-1],aura=min_bin_size)


            print "f",bin_bound_visualize(km_points,group[0],group[-1],aura=min_bin_size)
            print "g",bin_bound_visualize(expanded_group,group[0],group[-1])

            #calculate sum of variance from closest centroids
            dists = [[]]
            nearest = 0
            for cp in group:
                while nearest+1<len(km_points) and abs(cp-km_points[nearest+1]) < abs(cp-km_points[nearest]):
                    nearest+=1
                    dists.append([])
                dists[-1].append(abs(cp-km_points[nearest]))
            average_average_group_dist = sum(sum(ds)/float(len(ds)) for ds in dists)
            print "S",average_average_group_dist

            #add the solution to the list of possibilities
            solution_list.append((average_average_group_dist,km_points))

        #choose the best solution (that with the least variance)    
        best_score,best_bounds = sorted(solution_list)[0]
        #add the bounds determined for this group to the list of all determined bounds
        bin_bounds+=best_bounds

    #sort the determined bounds
    bin_bounds.sort()
    print "\nBin Bounds:", bin_bounds
    
    #Using the determined bounds, genotype each bin across RILs. Bins are genotyped as whichever genotype occupies the most 'area' in the bin.
    bin_genotypes = OrderedDict()
    for line in line_cps:
        bin_genotypes[line] = []

        line_cp_locs = list(enumerate((cp for i,cp in enumerate(line_cps[line]) if i%2 == 0)))
        scan_lower_i,scan_lower_bound = line_cp_locs[0]
        scan_upper_i,scan_upper_bound = line_cp_locs[0]

        bin_iter = enumerate(bin_bounds)
        next(bin_iter)
        for bin_i,bin_upper_bound in bin_iter:
            bin_lower_bound = bin_bounds[bin_i-1]

            while scan_lower_i+1<len(line_cp_locs) and line_cp_locs[scan_lower_i+1][1]<=bin_lower_bound:
                scan_lower_i,scan_lower_bound = line_cp_locs[scan_lower_i+1]
            while scan_upper_bound < bin_upper_bound and scan_upper_i+1<len(line_cp_locs):
                scan_upper_i,scan_upper_bound = line_cp_locs[scan_upper_i+1]

            scan_range_contents = line_cps[line][scan_lower_i*2:scan_upper_i*2+1]
            bin_contents = [bin_lower_bound]+scan_range_contents[1:-1]+[bin_upper_bound]

            genotype_weights = {}
            for seg_i in range(1,len(bin_contents),2):
                if not bin_contents[seg_i] in genotype_weights: genotype_weights[bin_contents[seg_i]] = 0
                genotype_weights[bin_contents[seg_i]] += bin_contents[seg_i+1]-bin_contents[seg_i-1]
            bin_genotypes[line].append(max(((genotype,genotype_weights[genotype]) for genotype in genotype_weights),key=lambda x:x[1])[0])

    #Find the center of each bin.
    bin_centers = [int((bin_bounds[i-1]+bin_bounds[i])/2) for i in range(1,len(bin_bounds))]

    #Remove any adjacent bins that are identical across all lines.
    for i in range(len(bin_centers)-1,0,-1):
        if all(bin_genotypes[line][i-1]==bin_genotypes[line][i] for line in bin_genotypes):
            print "Combined Identical Bins.", bin_bounds[i-1],">",bin_bounds[i],"<",bin_bounds[i+1]
            del bin_centers[i]
            del bin_bounds[i]
            bin_centers[i-1] = (bin_bounds[i-1]+bin_bounds[i])/2
            for line in bin_genotypes:
                del bin_genotypes[line][i]
    
    #Save the results!
    with open(output_path,"w") as outfile:
        if binmap_id: outfile.write("##binmap id,"+",".join([binmap_id]*len(bin_centers))+"\n")
        outfile.write("##bin start,"+",".join([str(i) for i in bin_bounds[:-1]])+"\n")
        outfile.write("##bin end,"+",".join([str(i) for i in bin_bounds[1:]])+"\n")
        outfile.write("bin center,"+",".join([str(i) for i in bin_centers])+"\n")
        for line in bin_genotypes:
            outfile.write("%s,%s\n" % (line,",".join(str(x) for x in bin_genotypes[line])))

def crosspoint_avg(weights,cp_list,chrom_len):
    '''Calculates the average for a list of crosspoints. Or, if one is at either end of the chromosome, return the near end of the chrom.'''
    if len(cp_list)<1:
        return float('nan')
    avg = sum(k*weights[k] for k in cp_list)/float(sum(weights[k] for k in cp_list))
    if math.isnan(avg) and (0 in cp_list or chrom_len in cp_list):
            avg = 0 if 0 in cp_list else chrom_len
    return int(avg)
                
def bin_bound_visualize(cps,begin,end,bins=75,aura=0):
    '''Returns a string of characters showing spacing, counts, and overlap bewteen bin bounds or crosspoints.'''
    empty_sym,aura_sym,aura_end_sym,aura_overlap_sym = (" ","-","|","x")
    bin_size = (end-begin)/float(bins-1)
    binned = [0]*(bins)
    aura_size = int(aura//bin_size//2 if aura!=0 else 0)
    for cp in cps:
        bin_loc = int((cp-begin)//bin_size)
        binned[bin_loc] += 1
    loc_map = [empty_sym]*(bins)
    for i in range(len(binned)):
        if binned[i]>0:
            if not binned[i]>10:
                loc_map[i] = str(binned[i])
            else:
                loc_map[i] = "G"
            if aura_size>0:
                right_end = i+aura_size
                left_end = i-aura_size
                if right_end<len(loc_map):
                    if loc_map[right_end] in (empty_sym,aura_end_sym): 
                        loc_map[right_end]=aura_end_sym
                    elif loc_map[right_end] in (aura_sym,aura_overlap_sym):
                        loc_map[right_end] = aura_overlap_sym
                if left_end>=0:
                    if loc_map[left_end] in (empty_sym,aura_end_sym): 
                        loc_map[left_end]=aura_end_sym
                    elif loc_map[left_end] in (aura_sym,aura_overlap_sym):
                        loc_map[left_end] = aura_overlap_sym
            for j in range(i-aura_size+1,i):
                if j>=0:
                    if loc_map[j]==empty_sym:
                        loc_map[j]=aura_sym
                    elif loc_map[j] in (aura_sym,aura_end_sym):
                        loc_map[j]=aura_overlap_sym
            for j in range(i+1,i+aura_size):
                if j<len(binned):
                    if loc_map[j]==empty_sym:
                        loc_map[j]=aura_sym
                    elif loc_map[j] in (aura_sym,aura_end_sym):
                        loc_map[j]=aura_overlap_sym
    return "".join(loc_map)