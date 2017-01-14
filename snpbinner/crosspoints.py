'''Uses genotyped SNP data to identify likely crossover points. (See README for more information)'''
from math import log

def parser(parser_add_func,name):
    '''Sets up an arguement parser for this module. Note the `dest` values match the arguements in the `run` function.'''
    p = parser_add_func(name,description="")
    #required
    p.add_argument("-i","--input",       metavar="PATH",  dest='input_file',                           required=True, help="Path to a SNP TSV.")
    p.add_argument("-o","--output",      metavar="PATH",  dest='output_file',                          required=True, help="Path for the output CSV.")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("-m","--min-length",  metavar="INT",   dest='min_state_length',      type=int,                     help="Minimum distance between crosspoints in basepairs. Cannot be used with `min-ratio`.")
    g.add_argument("-r","--min-ratio",   metavar="FLOAT", dest='min_state_ratio',       type=float,                   help="Minimum distance between crosspoints as a ratio. (0.01 would be 1%% of the chromosome.) Cannot be used with `min-length`.")
    #optional
    p.add_argument("-c","--cross-count", metavar="FLOAT", dest='predicted_cross_count', type=float,    default=4,     help="Used to calculate transition probability. The state transition probability is this value divided by the chromosome length. (default = 4)")
    p.add_argument("-l","--chrom-len",   metavar="INT",   dest='chrom_len',             type=int,      default=0,     help="The length of the chromosome/scaffold which the SNPs are on. If no length is provided, the last SNP is considered to be the last site on the chromosome.")
    p.add_argument("-p","--homogeneity", metavar="FLOAT", dest='predicted_homogeneity', type=float,    default=0.9,   help="Used to calculate emmision probabilities. For example, if 0.9 is used, it is predicted that a region genotyped as b would contain 90%% b-genotyped SNPs. (default = 0.9)")
    return p

def run(input_file,output_file,predicted_homogeneity,predicted_cross_count,chrom_len,min_state_length=None,min_state_ratio=None):
    '''Runs the module'''

    #get file statistics
    individual_count,snp_count,auto_chrom_len = get_file_stats(input_file)
    if chrom_len==0: 
        chrom_len=auto_chrom_len

    #determine the min distance between two crossover points
    if min_state_ratio!=None:
        min_state_length = chrom_len*float(min_state_ratio)

    # clear output file
    with open(output_file,"w") as outfile:
        pass

    # contruct the Hidden Markov Models (HMM) used to predict states
    crosspoint_probability = predicted_cross_count/float(chrom_len)
    intra_p = predicted_homogeneity # emmision probability of the genotype corresponding to the state
    inter_p = 1-intra_p                   # emmision probability of the opposite genotype

    # create a HMM that registers heterogeneous regions
    error = crosspoint_probability/2.0  # transition probability for other 2 states
    crrct = 1-(error*2)                 # transition probability for same state
    hmm_all = HMM(
        states      = ["a", "b", "h"],
        priors      = [0.3, 0.3, 0.3],
        transition = [[crrct, error, error],
                      [error, crrct, error],
                      [error, error, crrct]],
        observable  = ["a", "b", "h"],
        emission   = [[intra_p, inter_p, inter_p],   #Note that the probabilities dont sum to one in a&b emmision, as non-state SNPs are simply considered errors
                      [inter_p, intra_p, inter_p],
                      [0.5,     0.5,     inter_p*2]] #this configuration for h emmision has good performance, but does not have statistical meaning.
        )

    # create a HMM that does not register heterogeneous regions, this is used to split hetero regions more appropriately when it is required
    error = crosspoint_probability
    crrct = 1-error
    hmm_nohet = HMM(
        states      = ["a", "b"],
        priors      = [0.5, 0.5],
        transition = [[crrct, error],
                      [error, crrct]],
        observable  = ["a", "b"],
        emission   = [[intra_p, inter_p],
                      [inter_p, intra_p]] 
        )

    #Runs the crosspoint identification on the input file using the created HMMs and writes them to a file
    for i in range(0,individual_count):
        snplist,name = read_column(input_file,i)
        # get the crosspoints
        cp = find_crosspoints(
            snplist          = snplist,
            min_state_length = min_state_length,
            chrom_length     = chrom_len,
            hmm_nohet         = hmm_nohet,
            hmm_all          = hmm_all)
        with open(output_file,"a") as outfile:
            outfile.write(",".join([name]+[str(n) for n in cp])+",\n")
        print (name)

def find_crosspoints(snplist, min_state_length, chrom_length, hmm_nohet, hmm_all):
    '''Identifies crosspoints using two HMMs.'''

    # run HMM.gapped_viterbi to find state regions regardless of length
    cross_points  = hmm_all.gapped_viterbi(snplist)

    # run HMM.gapped_viterbi to find state regions regardless of length, ignoring the possibility of hetero regions.
    no_hetero_cross_points = hmm_nohet.gapped_viterbi(snplist)

    #Replace any heterogenous regions which are below the minimum state length with the most likely homogenous state regions
    for i in range(2,len(cross_points)-1,2):
        if cross_points[i]-cross_points[i-2]<min_state_length:
            if cross_points[i-1] == hmm_all.states[-1] and (hmm_all.states[-1] not in (cross_points[i-3],cross_points[i+1])):
                start,stop = None,None
                #find the region in the non-heterogeneous crosspoints that overlap with the too-short heterogeneous region
                for index,val in list(enumerate(no_hetero_cross_points))[0::2]:
                    if cross_points[i]<=val and stop==None:
                        stop = index+1
                        break
                    if cross_points[i-2]>=val:
                        start = index
                if not stop:
                    stop = no_hetero_cross_points[-1]
                #Replace the too-short heterogeneous region crosspoints with the identified non-heterogeneous region crosspoints
                cross_points[i-1:i] = no_hetero_cross_points[start+1:stop-1]

    #Identify all regions which are below the min state length and group contiguous segments of theses too-short regions.
    too_short = []
    for i in range(2,len(cross_points),2):
        if cross_points[i]-cross_points[i-2]<min_state_length:
            if len(too_short)>0 and too_short[-1][-1][1] == i-2:
                too_short[-1].append((i-2,i))
            else:
                too_short.append([(i-2,i)])

    # The algorithm used below requires two surrounding regions to function. However, It is both likely and possible that too-short regions will be at either end of the chromosome. In order to compensate for this, if an edge region is too short, it is marked as such but ignored in the algorithm. Instead, they will be assigned to their neighbors' genotypes after the algorithm runs.
    front_skipped,back_skipped = False,False
    if too_short and too_short[0][0][0]==0:
        front_skipped = True
        del too_short[0][0]
        if not too_short[0]:
            del too_short[0]
    if too_short and too_short[-1][-1][1]==len(cross_points)-1:
        back_skipped = True
        del too_short[-1][-1]
        if not too_short[-1]:
            del too_short[-1]


    # For each contiguous group of too-short regions, assign them genotypes. The algorithm is designed to minimize the number of mismatches between identified state and assigned state. In order to do this genotypes are assigned to contigous groups of too-short regions using the following rules, in order of priority where each rule is only applied when the conditions of those above it are not fullfilled:
    #   1. If a contiguous group of too-short regions is long enough to be its own acceptably-long region, it will be treated as such and assigned the most likely genotype using the 3-state HMM.
    #   2. If the first or last too-short region is neighboring to an acceptably-long region of the same genotype, it can be considered part of that region, and removed from the group.
    #   3. If the first or last too-short region is neighboring an acceptably-long heterogenous region, it will be assigned the heterogenous genotype and removed from the group as per Rule 2.
    #   4. If neither the first or last too-short region is neighboring a heterogenous or same-genotype region, the shortest of those two regions will be assigned to the same genotype as the acceptably-long region neighboring it and then removed as per Rule 2.
    #These steps repeat until each group is empty.
    for group in too_short:
        while len(group)>0:

            # Check if the leftmost (first) too-short region in a group matches the region to its left
            left_match = cross_points[group[0][0]-1] == cross_points[group[0][1]-1]

            # Check if the rightmost (last) too-short region in a group matches the region to its right
            right_match = cross_points[group[-1][1]+1] == cross_points[group[-1][1]-1] 

            if left_match or right_match: 
                # The first or last too-short region is neighboring an acceptably-long region of the same genotype
                if len(group)==1: 
                    # The only too-short region in the group already matches one of its neighbors, finished with group.
                    break
                if left_match:  
                    # The first too-short region matches its neighbor, so remove the first region from the group
                    group[:] = group[1:]
                if right_match: 
                    # The last too-short region matches its neighbor, so remove the last region from the group
                    group[:] = group[:-1]

            else: 
                # The genotypes do not match at either end of the group, so we need to change some genotypes!
                if cross_points[group[-1][1]]-cross_points[group[0][0]] >= min_state_length:
                    # The length of the whole group is enough to be its own region. So, treat all regions as one and assign the most likely genotype.
                    
                    #finds the SNPs represented by the region and runs them through the 3-state HMM (hmm_all)
                    first_snp = None
                    last_snp = None
                    i = 0
                    while (first_snp == None or last_snp == None) and i<len(snplist):
                        if snplist[i][0] == cross_points[group[0][0]]: first_snp = i
                        if snplist[i][0] == cross_points[group[-1][1]]: last_snp = i
                        i+=1
                    group_snps = snplist[first_snp:last_snp+1]
                    most_likely = hmm_all.get_most_likely_state(group_snps)
                    for reg in group:
                        cross_points[reg[1]-1] = most_likely
                    group = []
                    continue

                
                if hmm_all.states[-1] in (cross_points[group[0][0]-1],cross_points[group[-1][1]+1]):
                    # The first or last too-short region is neighboring an acceptably-long heterogenous region, it will be assigned the heterogenous genotype and removed from the group as per Rule 2.
                    if cross_points[group[0][0]-1] == hmm_all.states[-1]: 
                        # Hetero region on the left
                        cross_points[group[0][1]-1] = hmm_all.states[-1]
                        group[:] = group[1:]
                        if len(group)==0: continue
                    if cross_points[group[-1][1]+1] == hmm_all.states[-1]: 
                        # Hetero region on the right
                        cross_points[group[-1][1]-1] = hmm_all.states[-1]
                        group[:] = group[:-1]

                else: 
                    # There are no heterogenous matching surrounding regions, so assign the shortest of the two edge too-short regions to the genotype of its neighboring acceptable region.
                    left_size = cross_points[group[0][1]]-cross_points[group[0][0]]
                    right_size = cross_points[group[-1][1]]-cross_points[group[-1][0]]
                    if left_size<right_size:
                        cross_points[group[0][1]-1] = cross_points[group[0][0]-1]
                        group[:] = group[1:]
                    else:
                        cross_points[group[-1][1]-1] = cross_points[group[-1][1]+1]
                        group[:] = group[:-1]

    # Now, if a front or back too-short region was present, assign it to it's neighbor's genotype.
    if len(cross_points)>3:
        if front_skipped: cross_points[1] = cross_points[3]
        if back_skipped: cross_points[-2] = cross_points[-4]

    # Collapse the crosspoints so that there are none between regions of the same genotype.
    for i in range(len(cross_points)-3,1,-2):
        if cross_points[i-1] == cross_points[i+1]:
            del cross_points[i:i+2]

    # Check once more that the ends are not too short, if they are, collapse them.
    while len(cross_points)>2 and cross_points[2]<min_state_length:
        del cross_points[1:3]
    while len(cross_points)>2 and cross_points[-1]-cross_points[-3]<min_state_length:
        del cross_points[-3:-1]

    # Lastly, set the final length to the chrom_length
    cross_points[-1] = chrom_length

    return cross_points


class HMM(object):
    """A basic hidden markov model class for running the HMM.gapped_viterbi method"""
    def __init__(self,states,priors,transition,observable,emission):
        self.states = states
        self.priors = [log(n) for n in priors]
        self.transition = [[log(n) for n in line] for line in transition]
        self.observable = {pair[1]:pair[0] for pair in enumerate(observable)}
        self.emission = [[log(n) for n in line] for line in emission]

    def get_most_likely_state(self,all_obs_points):
        '''Calculates the probability a list of SNPs is of each single state then returns the most probable state.'''
        obs_points = [snp for snp in all_obs_points if snp[1] in self.observable]
        enum_states = list(enumerate(self.states))
        state_sums = ((sum(self.emission[i][self.observable[snp[1]]] for snp in obs_points),state) for i,state in enum_states)
        best_state = max(state_sums,key=lambda x: x[0])[1]
        return best_state
        
    def gapped_viterbi(self,all_obs_points):
        """A modified viterbi algorithm which outputs crosspoints
        Done using log(prob).

        obs_points is a sorted list of tuples conataing the position of 
        the observation in the first index and the observation in the second"""
        vtrb = []
        trace_graph = {}

        obs_points = [snp for snp in all_obs_points if snp[1] in self.observable]

        chrom_size = obs_points[-1][0]-obs_points[0][0]

        enum_states = list(enumerate(self.states))

        for t,snp in enumerate(obs_points):
            obs = snp[1]
            vtrb.append([])
            if t==0:
                for i,state in enum_states:
                    vtrb[0].append(self.priors[i] + self.emission[i][self.observable[obs]])
                    trace_graph[(0,i)] = None
            else:
                gap_size = snp[0] - obs_points[t-1][0]
                for i,state in enum_states:
                    gap_transition_adjustment = log(gap_size)
                    max_prob,max_prev = max((vtrb[-2][p] + (self.transition[p][i]+gap_transition_adjustment) + self.emission[i][self.observable[obs]], p) for p,prev_state in enum_states)
                    vtrb[t].append(max_prob)
                    trace_graph[(t,i)] = (t-1,max_prev)

        max_final,last = max((vtrb[-1][i],(len(vtrb)-1,i)) for i,state in enum_states)
        state_dets = []
        prev = last
        bp = []
        while prev!=None:
            state_dets.append(self.states[prev[1]])
            if trace_graph[prev] and trace_graph[prev][1]!=prev[1]:
                bp.append(prev[0])
            prev = trace_graph[prev]
        state_dets[:] = state_dets[::-1]
        bp[:] = bp[::-1]

        crosspoint_list = [0,state_dets[0]]
        for crosspoint in bp:
            crosspoint_list.append(obs_points[crosspoint][0])
            crosspoint_list.append(state_dets[crosspoint])
        crosspoint_list.append(obs_points[-1][0])

        return crosspoint_list

def read_column(filename, col, filter=True):
    '''Reads a single column with name header out of a TSV document '''
    with open(filename, "r") as f:
        snp_list = []
        name = ""
        title_line = ""
        while title_line.strip()=="" or title_line.startswith("#"):
            title_line = f.readline().split("#")[0] # get a new line and remove any end-of-line comments
        name = title_line.strip().split('\t')[col+2]
        for line in f:
            line = line.strip()
            if line!="" and not line.startswith("#"):
                line = line.split("#")[0] # remove any end-of-line comments
                items = line.split('\t')
                if filter==False or not items[col+2].lower()=="-": 
                    snp_list.append((int(float(items[1])), items[col+2].lower()))
        return snp_list,name

def get_file_stats(filename):
    '''Counts coulumns[-1 header](individual_count), rows[-1 header](snp_count), and returns the last row header (last_index) from a TSV.'''
    individual_count = 0
    with open(filename, "r") as f:
        title_line = ""
        while title_line.strip()=="" or title_line.startswith("#"):
            title_line = f.readline()
        indvs = title_line.strip().split('\t')[2:]
        individual_count = len(indvs)
    snps = read_column(filename,0,filter=False)[0]
    last_index = snps[-1][0]
    snp_count = len(snps)
    return individual_count,snp_count,last_index