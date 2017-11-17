from new_hmm import HMM
from pathlib import Path

def main():
    train(Path("./ch01.csv"),Path("./ch01.cps.csv"))

# def train(input_path,output_path):
def train(genotypes_path,crosspoints_path):
    '''Runs the module'''
    training_hmm = HMM()
    with  crosspoints_path.open() as crosspoints, genotypes_path.open() as genotypes:
        genotype_headers = genotypes.readline().split("\t")
        genotype_column = 0
        for cpline in crosspoints:
            # go to top and skip first line
            genotypes.seek(0)
            genotypes.readline()
            
            # read this line's cps
            cells = cpline.strip().strip(",").split(",")
            cpline_name = cells[0].strip()
            cpbounds = (int(n) for n in cells[3::2]) # crosspoint bounds
            cpgenos = iter(cells[2::2]) # genotypes
                        
            # find corresponding genotype column
            while genotype_headers[genotype_column].strip()!=cpline_name \
                    and genotype_column<len(genotype_headers):
                genotype_column+=1
                
            snps = snpGen(genotypes,genotype_column)
            
            training_hmm.train_begin()
            upper_bound,state = next(cpbounds),next(cpgenos)
            last_locus = None
            for loc,obs in snps:
                if last_locus is None: last_locus = loc-1
                while loc>upper_bound:
                    upper_bound,state = next(cpbounds),next(cpgenos)
                training_hmm.train(obs,state)
                last_locus = loc
            training_hmm.train_end()
            print(cpline_name)
            
        print("\n",str(training_hmm))
        new_hmm = HMM()
        new_hmm.load(training_hmm.dump())
        print("\n",str(new_hmm))
            
def snpGen(gt_file,genotype_column):
    for line in gt_file:
        cells = line.split("\t")
        pos = int(cells[1].strip())
        genotype = cells[genotype_column].strip()
        if genotype!="-":
            yield (pos,genotype)
    

# def _parser(parser_add_func,name):
#     '''Sets up an arguement parser for this module. Note the arguement names match those in the `bins` function.'''
#     p = parser_add_func(name,description=__doc__)
#     p.add_argument("-i","--input",        metavar="PATH",  dest='input_path',   required=True, nargs='*', help="Path to a binmap file (similar to the output of bins)")
#     p.add_argument("-o","--output",       metavar="PATH",  dest='output_path',  required=True,            help="Path for the JSON output")
#     return p

if __name__ == '__main__':
    main()
