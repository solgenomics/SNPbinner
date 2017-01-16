'''This script provides a way to visualize the results and input to the `bins` and `crosspoints` scripts. It can accept three filetypes (SNP input TSV, crosspoint script output CSV, and bin script output CSV). It then parses the files and groups the data by RIL, creating an image for each. In each colored row of the resulting images, regions are colored red, green, or blue, for genotype a, heterozygous, or genotype b, respecively. The binamp is represented in gray with adjacent bins alternating dark and light. The script can accept any combination or number of files for each of the different types.'''
import sys
import os
from collections import OrderedDict
import math
import numpy as np
from PIL import Image,ImageDraw, ImageFont

def parser(parser_add_func,name):
    p = parser_add_func(name,description="")
    p.add_argument("-o","--out",         metavar="PATH", dest="out_folder",      required=True,                 help="Folder to which the resulting images should be saved.")
    p.add_argument("-b","--bins",        metavar="PATH", dest="binned_path",     action='append', default = [], help="bins output file to be added to the visualization.")
    p.add_argument("-c","--crosspoints", metavar="PATH", dest="crosspoint_path", action='append', default = [], help="crosspoints output file to be added to the visualization.")
    p.add_argument("-s","--snps",        metavar="PATH", dest="snp_path",        action='append', default = [], help="SNP (crosspoints input file) file to be added to the visualization.")
    return p

def run(snp_path,crosspoint_path,binned_path,out_folder):
    snp_data, crosspoint_data, binned_data = OrderedDict(),OrderedDict(),OrderedDict()
    maximum_index = 0
    columns = 1000
    auto_set_columns = False
    min_columns = 500

    if snp_path:
        for file_path in snp_path:
            name = "SNPs: %s" % (file_path.strip().split("/")[-1])
            snp_data[name] = OrderedDict()
            with open(file_path) as snp_file:
                indv_names = snp_file.readline().strip().split("\t")[2:]
                for indv in indv_names:
                    snp_data[name][indv.strip()] = {"list":[]}
                for line in snp_file:
                    marker,position,gntyps = line.strip().split("\t",2)
                    for indv_i,genotype in enumerate(gntyps.split("\t")):
                        snp_data[name][indv_names[indv_i]]["list"].append([int(position),genotype])
                        maxi = snp_data[name][indv_names[indv_i]]["list"][-1][0]
                        if maxi>maximum_index:maximum_index=maxi
    if crosspoint_path:
        for file_path in crosspoint_path:
            name = "CPs: %s" % (file_path.strip().split("/")[-1])
            crosspoint_data[name] = {"indvs":OrderedDict(),"min_size":None}
            with open(file_path) as crosspoint_file:
                for line in crosspoint_file:
                    indv,cps = line.strip().strip(",").split(",",1)
                    indv = indv.strip()
                    crosspoints = cps.split(",")
                    crosspoints[::2] = (int(x) for x in crosspoints[::2])
                    indv_min_size = min(crosspoints[i]-crosspoints[i-2] for i in range(2,len(crosspoints),2))
                    crosspoint_data[name]["indvs"][indv] = {"list":crosspoints}
                    if crosspoints[-1]>maximum_index:maximum_index=crosspoints[-1]
    if binned_path:
        for file_path in binned_path:
            name = "BINs: %s" % (file_path.strip().split("/")[-1])
            binned_data[name] = {"indvs":OrderedDict(),"bin_bounds":None,"min_bin_size":None}
            with open(file_path) as bin_file:
                line = bin_file.readline()
                while not line.startswith("##bin start"):
                    line = bin_file.readline()
                bin_starts = [int(i) for i in line.strip().split(',')[1:]]
                bin_ends = [int(i) for i in bin_file.readline().strip().split(',')[1:]]
                bin_centers = [int(i) for i in bin_file.readline().strip().split(',')[1:]]
                maxi = bin_ends[-1]
                if maxi>maximum_index:maximum_index=maxi
                binned_data[name]["bin_bounds"] = [0]+bin_ends
                binned_data[name]["min_bin_size"] = min(binned_data[name]["bin_bounds"][i]-binned_data[name]["bin_bounds"][i-1] for i in range(1,len(binned_data[name]["bin_bounds"])))
                for line in bin_file:
                    indv,gntyps = line.strip().split(",",1)
                    indv = indv.strip()
                    genotypes = gntyps.split(",")
                    all_bins = [None]*(len(genotypes)+len(binned_data[name]["bin_bounds"]))
                    all_bins[::2] = binned_data[name]["bin_bounds"]
                    all_bins[1::2] = genotypes
                    only_cp = all_bins[:2]
                    for i in range(3,len(all_bins),2):
                        if all_bins[i]!=only_cp[-1]:
                            only_cp+=all_bins[i-1:i+1]
                    only_cp.append(all_bins[-1])
                    binned_data[name]["indvs"][indv] = {"list":only_cp}
    column_size = (maximum_index)/float(columns-1)
    if auto_set_columns:
        if len(binned_data)>0:
            min_size = min(binned_data[name]["min_bin_size"] for name in binned_data)
            columns = maximum_index//min_size + 1
            column_size = maximum_index/float(columns-1)
        elif len(crosspoint_data)>0:
            min_size = min(binned_data[name]["min_size"] for name in binned_data)
            columns = maximum_index//min_size + 1
            column_size = maximum_index/float(columns-1)

    for name in snp_data:
        for indv in snp_data[name]:
            snp_data[name][indv]["dict"] = {}
            for snp in snp_data[name][indv]["list"]:
                col = int(math.floor(snp[0]/column_size))
                if snp[1] in genotype_values:
                    if not col in snp_data[name][indv]["dict"]:
                        snp_data[name][indv]["dict"][col] = []
                    snp_data[name][indv]["dict"][col].append(snp[1])
            snp_data[name][indv]["col_vals"] = [None]*columns
            for col in snp_data[name][indv]["dict"]:
                snp_data[name][indv]["col_vals"][col] = sum(genotype_values[x] for x in snp_data[name][indv]["dict"][col])/float(len(snp_data[name][indv]["dict"][col]))
    for cp_data_set in (crosspoint_data,binned_data):
        for name in cp_data_set:
            for indv in cp_data_set[name]["indvs"]:
                cp_data_set[name]["indvs"][indv]["list"][::2] = (int(math.floor(x/column_size)) for x in cp_data_set[name]["indvs"][indv]["list"][::2])
                collapsed_cp_cols = [cp_data_set[name]["indvs"][indv]["list"][0]]
                genotype_list = []
                for i in range(2,len(cp_data_set[name]["indvs"][indv]["list"]),2):
                    genotype_list.append(cp_data_set[name]["indvs"][indv]["list"][i-1])
                    if cp_data_set[name]["indvs"][indv]["list"][i]!=cp_data_set[name]["indvs"][indv]["list"][i-2]:
                        collapsed_cp_cols.append(sum(genotype_values[x] for x in genotype_list)/float(len(genotype_list)))
                        collapsed_cp_cols.append(cp_data_set[name]["indvs"][indv]["list"][i])
                        genotype_list = []
                if len(genotype_list)>0:
                    collapsed_cp_cols.append(sum(genotype_values[x] for x in genotype_list)/float(len(genotype_list)))
                    collapsed_cp_cols.append(cp_data_set[name]["indvs"][indv]["list"][-1])
                cp_data_set[name]["indvs"][indv]["col_vals"] = [None]*columns
                for i in range(2,len(collapsed_cp_cols),2):
                    cp_data_set[name]["indvs"][indv]["col_vals"][collapsed_cp_cols[i-2]:collapsed_cp_cols[i]] = (collapsed_cp_cols[i-1] for n in range(collapsed_cp_cols[i-2],collapsed_cp_cols[i]))

    for name in binned_data:
        binned_data[name]["bin_bounds"][:] = (int(math.floor(x/column_size)) for x in binned_data[name]["bin_bounds"][:])
        binned_data[name]["bin_col_vals"] = [None]*columns
        color_alt = 1
        for i in range(1,len(binned_data[name]["bin_bounds"])):
            binned_data[name]["bin_col_vals"][binned_data[name]["bin_bounds"][i-1]:binned_data[name]["bin_bounds"][i]] = [3+color_alt]*(binned_data[name]["bin_bounds"][i]-binned_data[name]["bin_bounds"][i-1])
            color_alt *= -1
    individuals = OrderedDict()
    for name in snp_data:
        for indv in snp_data[name]:
            if not indv in individuals: individuals[indv] = {"SNP_rows":{},"CP_rows":{},"BIN_rows":{},"BIN_maps":{}}
            individuals[indv]["SNP_rows"][name] = snp_data[name][indv]["col_vals"]
    for cp_data_set,cp_row_type in ((crosspoint_data,"CP_rows"),(binned_data,"BIN_rows")):
        for name in cp_data_set:
            for indv in cp_data_set[name]["indvs"]:
                if not indv in individuals: individuals[indv] = {"SNP_rows":{},"CP_rows":{},"BIN_rows":{},"BIN_maps":{}}
                individuals[indv][cp_row_type][name] = cp_data_set[name]["indvs"][indv]["col_vals"]
    for name in binned_data:
        for indv in individuals:
            if not indv in individuals: individuals[indv] = {"SNP_rows":{},"CP_rows":{},"BIN_rows":{},"BIN_maps":{}}
            individuals[indv]["BIN_maps"][name+" (binmap)"] =binned_data[name]["bin_col_vals"]
    image_builders = []
    row_order = [item.split("/")[-1] for item in sys.argv[2:] if not item.startswith("-")]
    for indv in individuals:
        row_list = []
        if not indv in individuals: individuals[indv] = {"SNP_rows":{},"CP_rows":{},"BIN_rows":{},"BIN_maps":{}}
        for name in individuals[indv]["SNP_rows"]:
            row_list.append((name,individuals[indv]["SNP_rows"][name]))
        for name in individuals[indv]["CP_rows"]:
            row_list.append((name,individuals[indv]["CP_rows"][name]))
        for name in individuals[indv]["BIN_rows"]:
            row_list.append((name,individuals[indv]["BIN_rows"][name]))
        for name in individuals[indv]["BIN_maps"]:
            row_list.append((name,individuals[indv]["BIN_maps"][name]))
        row_list.sort(key=lambda row:next((i for i, x in enumerate(row_order) if  row[0].split(":",1)[-1].strip().startswith(x)), -1))

        builder = Indvidual_image_builder(indv)
        for row in row_list:
            builder.add_row(*row)
        image_builders.append(builder)

    out_folder = out_folder.rstrip("/")+"/   "
    for builder in image_builders:
        builder.save(out_folder+builder.name+".png")
        




def rgba_hex_val(r,g,b,a=255): return reduce(lambda x,y:(x<<8)|y,(a,b,g,r),0)

genotype_values = {"a":1,"h":0,"b":-1}
genotype_colors = {"a":[255,0,0],"b":[0,0,255],"h":[0,255,0],"no_data":[255,255,255]}
bin_map_values = {"bmc1":2,"bmc2":4}
bin_map_colors = {"bmc1":[107,107,107],"bmc2":[147,147,147]}

def val_to_color(val):
    if val!=None and val<=1 and val>=0:
        ratio = val
        return rgba_hex_val(*[int(i*ratio+j*(1-ratio)) for i,j in zip(genotype_colors["a"],genotype_colors["h"])])
    elif val!=None and val>=-1 and val<0:
        ratio = abs(val)
        return rgba_hex_val(*[int(i*ratio+j*(1-ratio)) for i,j in zip(genotype_colors["b"],genotype_colors["h"])])
    elif val==bin_map_values["bmc1"]:
        return rgba_hex_val(*bin_map_colors["bmc1"])
    elif val==bin_map_values["bmc2"]:
        return rgba_hex_val(*bin_map_colors["bmc2"])
    return rgba_hex_val(*(genotype_colors["no_data"]))

class Indvidual_image_builder(object):
    title_height = 25
    padding = 1
    def __init__(self, name):
        self.name = name.strip()
        self.length = None
        self.row_list = OrderedDict()
    def add_row(self,row_name,row_col_cols):
        self.row_list[row_name] = [val_to_color(x) for x in row_col_cols]
        if self.length==None: self.length = len(self.row_list[row_name])
    def __str__(self):
        strin = self.name+"\n"
        for row_name in self.row_list:
            strin+="(%s)%s\n" % (row_name,self.row_list[row_name])
        return strin
    def save(self,file_name):
        sys.stderr.write("Saving %s..." % file_name)
        section_height = int(1/10.0 * self.length)
        w = self.length
        h = len(self.row_list)*section_height+self.title_height
        img = np.empty((w,h),np.uint32)
        img.shape=h,w
        img[0:h,0:w] = 0xffffffff
        for rn_i,row_name in enumerate(self.row_list):
            row_top = rn_i*section_height
            row_bot = row_top+section_height-self.padding
            for line in img[row_top:row_bot]:
                line[:] = self.row_list[row_name]
        pil_im = Image.frombuffer('RGBA',(w,h),img,'raw','RGBA',0,1)
        draw = ImageDraw.Draw(pil_im)
        font = ImageFont.truetype(os.path.join(os.path.dirname(__file__), "package_data/Tuffy.ttf"),size=14)
        for rn_i,row_name in enumerate(self.row_list):
            binned_label_middle = (rn_i+0.5)*section_height
            lw,lh = font.getsize(row_name)
            binned_label_top = binned_label_middle - lh//2
            binned_label_bottom = binned_label_top + lh
            binned_label_left = self.length//2 - lw//2
            binned_label_right = binned_label_left + lw
            text_loc = (self.length//2-lw//2,(binned_label_middle-lh//2))
            box_bounds = ((binned_label_left-3,binned_label_top),(binned_label_right+2,binned_label_bottom))
            draw.rectangle(box_bounds, fill=(255,255,255))
            draw.text(text_loc, row_name, (0,0,0),font=font)

        tw,th = font.getsize(self.name)
        draw.rectangle(((0,h-self.title_height),(w,h)), fill=(0,0,0))
        draw.text((self.length//2-tw//2,(h-(self.title_height//2)-th//2)), self.name, (255,255,255),font=font)

        pil_im.save(file_name,optimize=True)
        sys.stderr.write(" Done.\n")