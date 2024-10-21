MAIN CLASS (FI-DE)
import math
import numpy as np
import random
import pandas as pd
from statistics import mean

def DE(biclusters, biclu, best):
    dimen = len(biclusters)
    F = 0.2
    CR = 0.5
    abc = random.sample(biclusters,2)
    q = random.randint(1,dimen)
    for k in range(dimen):
        RK = random.uniform(0,1)
        if RK<CR or k==q:
            mutate = np.array(best) + F*(np.array(abc[0])-np.array(abc[1]))
            mutate_round = roundoff(mutate).tolist()
            new_bic = [round(x) for x in mutate_round]
        else:
            new_bic = biclu
    return new_bic

def hamming(biclu_1, biclu_2):
    diff = 0
    for ch1, ch2 in zip(biclu_1, biclu_2):
        if ch1 != ch2:
            diff += 1
    return diff    
def volume(biclusters,n):
    volumes = []
    for line in biclusters:
        vol = line[:n].count(1)*line[n:].count(1)
        volumes.append(vol)
    return volumes
def roundoff(bics):
    for j, item in enumerate(bics):
        if item < 0.95:
            bics[j] = 0
        else:
            bics[j] = 1
    return bics
def MSR(bic):
    bicluster_residue = []
    aiJ_list = []
    aIj_list = []

    #calculate the mean of each row in the bicluster
    for row in bic:
        row_mean = mean(row)
        aiJ_list.append(row_mean)

    #calculate the mean of each column in the bicluster
    bicluster_transpose = np.array(bic).T.tolist()
    for column in bicluster_transpose:
        column_mean = mean(column)
        aIj_list.append(column_mean)
    
    aIJ = mean(aiJ_list)    #mean of the bicluster

    #calculate the residue of each element in the bicluster
    for i, row in enumerate(bic):
        row_residue = []
        for j, item in enumerate(row):
            res = item + aIJ - aiJ_list[i] - aIj_list[j]
            row_residue.append(res)
        bicluster_residue.append(row_residue)

    #calculate the MSR of the bicluster
    residue_sum = 0.0
    I = len(bicluster_residue)
    J = len(bicluster_residue[0])

    for row in bicluster_residue:
        row_sum = 0.0
        for item in row:
            row_sum = row_sum + item**2
        residue_sum = residue_sum + row_sum
    msr = residue_sum/(I*J)
    
    #calculate the variance of the bicluster
    variance_sum = 0.0
    for i, row in enumerate(bic):
        sum_vari = 0.0
        for item in row:
            sum_vari += (item - aiJ_list[i])**2
        row_variance = sum_vari/J
        variance_sum += row_variance
    bi_variance = variance_sum/I

    #calculate the bicluster fitness
    if bi_variance == 0:
        bi_fitness = msr
    else:
        bi_fitness = msr + 1/bi_variance
    return msr, bi_variance, bi_fitness

def extract_bicluster(input_data, new_bic):
    gindex = [i for i, item in enumerate(new_bic[:15380]) if item==1]
    cindex = [i for i, item in enumerate(new_bic[15380:]) if item==1]
    bic = []
    for item in gindex:
        row = []
        for line in cindex:
            row.append(input_data[item][line])
        bic.append(row)
    return bic
#exclude the same biclusters using MSR
def exclude_same_biclusters(MSRs, var, fit, biclusters):
    for i, item1 in enumerate(MSRs):
        for j, item2 in enumerate(MSRs):
            if item1 == item2:
                if biclusters[i] == biclusters[j]:
                    del(biclusters[j])
                    del(MSRs[j])
                    del(var[j])
                    del(fit[j])
    return MSRs, biclusters       
def sort_biclusters(MSRs, var, fit, biclusters):
    zipped_pairs = zip(MSRs, var, fit, biclusters)
    sorted_pairs = sorted(zipped_pairs)
    tuple_pairs = zip(*sorted_pairs)
    sorted_msr, sorted_var, sorted_fit, sorted_bic = [list(tuple) for tuple in  tuple_pairs]
    return sorted_msr, sorted_var, sorted_fit, sorted_bic
def annd(f,bic1):
    dim = len(bic1)
    and_result = []
    for i in range(dim):
        if bic1[i] == 1 and f == 1:
            and_result.append(1)
        else:
            and_result.append(0)   
    return and_result
def xor(bic1,bic2):
    dim = len(bic1)
    xor_result = []
    for i in range(dim):
        if bic1[i]==bic2[i]:
            xor_result.append(0)
        else:
            xor_result.append(1)
    return xor_result
def union(bic1,bic2):
    dim = len(bic1)
    union_result = []
    for j in range(dim):
        if bic1[j] == 0 and bic2[j] == 0:
            union_result.append(0)
        else:
            union_result.append(1)   
    return union_result
def index(biclus):
    biclusters_gene_cond_count = []
    for item in biclus:
        gene_cond_count = []
        gene_count = item[:15380].count(1)
        cond_count = item[15380:].count(1)
        gene_cond_count.append(gene_count)
        gene_cond_count.append(cond_count)
        biclusters_gene_cond_count.append(gene_cond_count)
    return biclusters_gene_cond_count
#reading encoded bicluster file and normalized input data
df = pd.read_csv(r'C:/Users/THRIVESBIO/Documents/PD/Mrs Osuntokun/project/Biclustering/output/encoded_schisto.csv', header=None)
df2 = pd.read_csv(r'C:/Users/THRIVESBIO/Documents/PD/Mrs Osuntokun/project/Biclustering/output/norm_input_file_schisto.csv', header=None)
df.dropna(axis=0, how='all',inplace=True)
df2.dropna(axis=0, how='all',inplace=True)
biclust = df.values.tolist()
input_data = df2.values.tolist()

biclusters = []
for line in biclust:
    biclusters.append([int(x) for x in line])

#calculate MSR, gene variance and fitness of all seeds
MSRs = []
var = []
fit = []
for item in biclusters:
    seed_msr = extract_bicluster(input_data, item)
    MSRs.append(MSR(seed_msr)[0])
    var.append(MSR(seed_msr)[1])
    fit.append(MSR(seed_msr)[2])

print(len(MSRs))
MSRs, biclusters = exclude_same_biclusters(MSRs, var, fit, biclusters)
print()
print(len(MSRs))

MSRs, var, fit, biclusters = sort_biclusters(MSRs, var, fit, biclusters)
biclusters_dimen = len(biclusters)
print(index(biclusters))
print(MSRs)
beta = 1.0
gamma = 5.0
alpha = 0.5
F = 1
#abc = random.sample(biclusters,2)
best_index = MSRs.index(min(MSRs))
best = biclusters[best_index]
for l in range(10):
    print('Epoch: ',l)
    for i in range(biclusters_dimen):
        ab = random.sample(biclusters,2)
        '''if i == biclusters_dimen-1:
            xor_res = xor(biclusters[0],biclusters[i])
        else:'''
        xor_res = xor(ab[0],ab[1])
        temp_bic = union(best, xor_res)
        new_bic = union(biclusters[i],temp_bic)
        new_bic_data = extract_bicluster(input_data, new_bic)
        new_bic_MSR = MSR(new_bic_data)[0]
        new_bic_var = MSR(new_bic_data)[1]
        new_bic_fit = MSR(new_bic_data)[2]
        biclusters[i] = new_bic
        MSRs[i] = new_bic_MSR
        var[i] = new_bic_var
        fit[i] = new_bic_fit
        MSRs, var, fit, biclusters = sort_biclusters(MSRs, var, fit, biclusters)
        ''' for j in range(biclusters_dimen):
            if MSRs[i] < MSRs[j]:
                #vol = biclusters[i].count(1)
                dist = hamming(biclusters[i],biclusters[j])
                attract = beta*math.exp(-gamma*dist**2)
                move = np.array(biclusters[i]) + attract*(np.array(biclusters[j]) - np.array(biclusters[i])) + alpha*(random.random() - 0.5)
                xor_res = xor(biclusters[j],biclusters[i])
                #and_res = annd(F, xor_res)
                new_bic = union(best, xor_res)
                #new_bic = union(biclusters[i], mutate)
                new_bic_data = extract_bicluster(input_data, new_bic)
                geneindex = [i for i, item in enumerate(new_bic[:721]) if item==1]
                condindex = [i for i, item in enumerate(new_bic[721:]) if item==1]
                if 1 not in geneindex or 1 not in condindex:
                   continue
                new_bic_MSR = MSR(new_bic_data)
                #new_vol = new_bic.count(1)
                #if new_bic_MSR < 300:
                biclusters[i] = new_bic
                MSRs[i] = new_bic_MSR
                 
            else:
               new_bic = DE(biclusters,biclusters[i],best)
               geneindex = [i for i, item in enumerate(new_bic[:721]) if item==1]
               condindex = [i for i, item in enumerate(new_bic[721:]) if item==1]
               if 1 not in geneindex or 1 not in condindex:
                   continue
               new_bic_data = extract_bicluster(input_data, new_bic)
               new_bic_MSR = MSR(new_bic_data)
               if new_bic_MSR > 300:
                   continue
               else:
                   biclusters[i] = new_bic
                   MSRs[i] =new_bic_MSR'''
    print(index(biclusters))            
    print(MSRs)
print(var[:15])
    #print(fit)
biclusters_index = []
for line in biclusters:
    biclu_index = []
    biclu_index.append([x for x,item in enumerate(line[:15380]) if item==1])
    biclu_index.append([x for x,item in enumerate(line[15380:]) if item==1])
    biclusters_index.append(biclu_index)
print(biclusters_index[:15])

Extract_gene names

import csv

path = r'C:/Users/THRIVESBIO/Documents/PD/Mrs Osuntokun/project/Biclustering/data/schisto_genenames.csv'

with open(path) as f:
    gene_names = csv.reader(f, delimiter = '\t')
    gnames = [line for sublist in gene_names for line in sublist]

path1 = r'C:/Users/THRIVESBIO/Documents/PD/Mrs Osuntokun/project/Biclustering/output/bicluster_gene_index_schisto.txt'
input_file  = open(path1)
bicluster_index = []
for line in input_file:
    gene_line = line.split(',')
    bicluster_index.append([int(i) for i in gene_line])

bicluster_genes = []
for line in bicluster_index:
    each_bicluster_genes = []
    for item in line:
        each_bicluster_genes.append(gnames[item])
    bicluster_genes.append(each_bicluster_genes)


for line in bicluster_genes:
    print(line)
    print()

Clustering. Py
import csv
import math
import numpy as np
from nltk.cluster import KMeansClusterer, cosine_distance
from statistics import mean

#function to recluster bicluster seed that has # of genes>15
def recluster_gene(pos,gene_data):
    biclut_array = [gene_data[x] for x in pos]
    gene_clusterer = KMeansClusterer(2, cosine_distance)
    gene_clusters = gene_clusterer.cluster(biclut_array, True, trace = True)
    group1 = []
    group2 =[]
    for i, item in enumerate(gene_clusters):
        if item ==0:
            group1.append(pos[i])
        else:
            group2.append(pos[i])
    return group1,group2

#function to recluster bicluster seed that has # of conditions>5
def recluster_condition(pos,cond_data):
    biclut_array = [cond_data[x] for x in pos]
    condition_clusterer = KMeansClusterer(2, cosine_distance)
    condition_clusters = condition_clusterer.cluster(biclut_array, True, trace = True)
    group1 = []
    group2 =[]
    for i, item in enumerate(condition_clusters):
        if item ==0:
            group1.append(pos[i])
        else:
            group2.append(pos[i])
    return group1,group2
#open input file
in_file = []
path = r'C:/Users/THRIVESBIO/Documents/PD/Mrs Osuntokun/project/Biclustering/data/schisto.csv'
with open(path) as f:
    input_file = csv.reader(f, delimiter = '\t')
    for line in input_file:
        gene_line = line[0].split(",")
        in_file.append([float(i) for i in gene_line])

#transformation of raw expression values with log function
normalized_file = []
for line in in_file:
	ll=[]
	for item in line:
            if item <=0:
                norm = 0
                ll.append(norm)
            else:
                norm = 100*(math.log(10**5*item))
                ll.append(norm)
	normalized_file.append(ll)
'''input_file  = open(path)
for line in input_file:
    gene_line = line.split()
    in_file.append([float(i) for i in gene_line])'''
    
#write normalized raw data
with open("C:/Users/THRIVESBIO/Documents/PD/Mrs Osuntokun/project/Biclustering/output/norm_input_file.csv","w") as f:
    wr = csv.writer(f)
    wr.writerows(in_file)  
    
#convert list data to array data
gene_normalized_array = np.array([np.array(i) for i in in_file])
condition_normalized_array = np.array([np.array(i) for i in in_file]).T

#cluster the data using gene, k=54(yeast), k=124(schisto)
gene_clusterer = KMeansClusterer(124, cosine_distance, repeats = 10)
gene_clusters = gene_clusterer.cluster(gene_normalized_array, True, trace = True)

#extract the genes positions for each cluster
gene_positions = []
for k in range(124):
    gene_pos = []
    for i, item in enumerate(gene_clusters):
        if item != k:
            continue
        else:
            gene_pos.append(i)
    gene_positions.append(gene_pos)

#print(gene_positions)
print(len(gene_positions))

#recluster seeds with higher number of genes
countg = 0
while countg<len(gene_positions):
    if len(gene_positions[countg])>80:
        grp1,grp2 = recluster_condition(gene_positions[countg],gene_normalized_array)
        del(gene_positions[countg])
        gene_positions.append(grp1)
        gene_positions.append(grp2)
    else:
        countg +=1
#print(gene_positions)
print(len(gene_positions))
#cluster the data using condition, k=4(yeast), k=5(shisto)
condition_clusterer = KMeansClusterer(5, cosine_distance, repeats = 5)
condition_clusters = condition_clusterer.cluster(condition_normalized_array, True, trace = True)

#extract the condition positions for each cluster
condition_positions = []
for k in range(5):
    condition_pos = []
    for i, item in enumerate(condition_clusters):
        if item != k:
            continue
        else:
            condition_pos.append(i)
    condition_positions.append(condition_pos)
print(condition_positions)  
print(len(condition_positions))

#recluster seeds with higher number of conditions
countc = 0
while countc<len(condition_positions):
    if len(condition_positions[countc])>5:
        grp1,grp2 = recluster_condition(condition_positions[countc],condition_normalized_array)
        del(condition_positions[countc])
        condition_positions.append(grp1)
        condition_positions.append(grp2)
    else:
        countc +=1

print(condition_positions)
print(len(condition_positions))

#extract biclusters
biclusters = []
biclusters_indices = []
for line in gene_positions:
    gene_index = line
    for row in condition_positions:
        each_bicluster = []
        each_bicluster_index = []
        for item in gene_index:
            each_row = []
            for column_index in row:
                each_row.append(gene_normalized_array[item,column_index])
            each_bicluster.append(each_row)
            each_bicluster_index.append(gene_index)
            each_bicluster_index.append(row)
        biclusters.append(each_bicluster)
        biclusters_indices.append(each_bicluster_index)
print(len(biclusters))

#MSR function
biclusters_MSR = []
for k, line in enumerate(biclusters):
    bicluster_residue = []
    aiJ_list = []
    aIj_list = []

    #calculate the mean of each row in the bicluster
    for row in line:
        row_mean = mean(row)
        aiJ_list.append(row_mean)

    #calculate the mean of each column in the bicluster
    bicluster_transpose = np.array(line).T.tolist()
    for column in bicluster_transpose:
        column_mean = mean(column)
        aIj_list.append(column_mean)

    aIJ = mean(aiJ_list)    #mean of the bicluster

    #calculate the residue of each element in the bicluster
    for i, row in enumerate(line):
        row_residue = []
        for j, item in enumerate(row):
            res = item + aIJ - aiJ_list[i] - aIj_list[j]
            row_residue.append(res)
        bicluster_residue.append(row_residue)

    #calculate the MSR of the bicuster
    residue_sum = 0.0
    I = len(bicluster_residue)
    J = len(bicluster_residue[0])

    for row in bicluster_residue:
        row_sum = 0.0
        for item in row:
            row_sum = row_sum + item**2
        residue_sum = residue_sum + row_sum
    MSR = residue_sum/(I*J)
    
    biclusters_MSR.append(MSR)
print(len(biclusters_MSR))
print(sorted(biclusters_MSR))
#remove biclusters with MSR>threshold(yeast=300, schisto=2100, lymphoma=1200)
new_biclusters = []
new_biclusters_MSR = []
new_biclusters_indices = []
for j, item in enumerate(biclusters_MSR):
    if biclusters_MSR[j]>=2100:
        continue
    else:
        new_biclusters.append(biclusters[j])
        new_biclusters_indices.append(biclusters_indices[j])
        new_biclusters_MSR.append(biclusters_MSR[j])
print()
print(len(new_biclusters))
print(len(new_biclusters_MSR))
#ENCODING
data_dimen = gene_normalized_array.shape
encoded_biclusters = []
result_file = open("C:/Users/THRIVESBIO/Documents/PD/Mrs Osuntokun/project/Biclustering/output/encoded.csv",'w')
wr = csv.writer(result_file, dialect='excel')
for line in new_biclusters_indices:
    g_index = line[0]
    c_index = line[1]
    gene_code = np.zeros(data_dimen[0], dtype = int)
    condition_code = np.zeros(data_dimen[1], dtype = int)
    for item in g_index:
        gene_code[item] = 1  
    for item in c_index:
        condition_code[item] = 1
    encoded_bi = gene_code.tolist() + condition_code.tolist()
    encoded_biclusters.append(encoded_bi)
    wr.writerow(encoded_bi)
print(sorted(new_biclusters_MSR))
print('Encoded file and biclusters written')


Biclustering. PY
import csv
import math
import numpy as np
from nltk.cluster import KMeansClusterer, cosine_distance
from statistics import mean


class file_input(object):
    
    def __init__(self,file_path):
        self.file_path=file_path
        
    #function that extracts the gene expression value from text input file
    def file_input1(self):
        path  = self.file_path
        in_file = []
        input_file  = open(path)
        for line in input_file:
            gene_line = line.split()
            in_file.append(gene_line)
        numeric_file = []
        for line in in_file:
            new = [float(i) for i in line]
            numeric_file.append(new)
        normalized_array = np.array([np.array(i) for i in numeric_file])   
        return normalized_array
    
    #function that extracts the gene expression value from csv input file
    def file_input2(self):
        path  = self.file_path
        in_file = []
        with open(path) as f:
            input_file = csv.reader(f, delimiter = '\t')
            for line in input_file:
                gene_line = line[0].split(",")
                in_file.append(gene_line)

        numeric_file = []
        for line in in_file:
            new = [float(i) for i in line]
            numeric_file.append(new)

        normalized_file = []
        for line in numeric_file:
            ll=[]
            for item in line:
                if item <=0:
                    continue
                norm = 100*(math.log(10**5*item))
                ll.append(norm)
            normalized_file.append(ll)
        normalized_array = np.array([np.array(i) for i in normalized_file])
        return normalized_array
    
    #function that performs the initial gene seed generation using kmeans
    def gene_seed(self, normalized_array, n):
        gene_seed_file = infile
        no_clusters = n
        clusterer = KMeansClusterer(no_clusters, cosine_distance)
        seed_gene = clusterer.cluster(gene_seed_file, True, trace = True)
        return seed_gene
        
    #function that performs the initial condition seed generation using kmeans
    def condition_seed(self, normalized_array, n):
        condition_seed_file = np.array(normalized_array).T
        no_clusters = n
        clusterer = KMeansClusterer(no_clusters, cosine_distance)
        seed_condition = clusterer.cluster(condition_seed_file, True, trace = True)
        return seed_condition

    #function that extracts the genes positions for each cluster
    def gene_positions(self, gene_clusters, n):
        gene_positions = []
        for k in range(n):
            gene_pos = []
            for i, item in enumerate(gene_clusters):
                if item != k:
                    continue
                else:
                    gene_pos.append(i)
            gene_positions.append(gene_pos)
        return gene_positions
    
    #function that reclusters seeds with higher number of genes
    def recluster_gene(gene_positions,gene_data):
        countg = 0
        while countg<len(gene_positions):
            if len(gene_positions[countg])>15:
                gene_biclut_array = [gene_data[x] for x in gene_positions[countg]]
                gene_clusterer = KMeansClusterer(2, cosine_distance)
                gene_clusters = gene_clusterer.cluster(gene_biclut_array, True, trace = True)
                group1 = []
                group2 = []
                for i, item in enumerate(gene_clusters):
                    if item ==0:
                        group1.append(gene_positions[countg][i])
                    else:
                        group2.append(gene_positions[countg][i])
                del(gene_positions[countg])
                gene_positions.append(group1)
                gene_positions.append(group2)
            else:
                countg +=1
    #function that extracts the genes positions for each cluster
    def condition_positions(self, condition_clusters, n):
        condition_positions = []
        for k in range(n):
            condition_pos = []
            for i, item in enumerate(condition_clusters):
                if item != k:
                    continue
                else:
                    condition_pos.append(i)
            condition_positions.append(condition_pos)
        return condition_positions
    
    #function that reclusters seeds with higher number of conditions
    def recluster_condition(condition_positions,cond_data):
        countc = 0
        while countc<len(condition_positions):
            if len(condition_positions[countc])>5:
                cond_biclut_array = [cond_data[x] for x in condition_positions[countc]]
                condition_clusterer = KMeansClusterer(2, cosine_distance)
                condition_clusters = condition_clusterer.cluster(cond_biclut_array, True, trace = True)
                group1 = []
                group2 = []
                for i, item in enumerate(condition_clusters):
                    if item ==0:
                        group1.append(condition_positions[countc][i])
                    else:
                        group2.append(condition_positions[countc][i])
                del(condition_positions[countc])
                condition_positions.append(group1)
                condition_positions.append(group2)
            else:
                countc +=1
    
    #function that extracts biclusters
    def extract_biclusters(self, gene_positions, condition_positions, normalized_array):
        biclusters = []
        biclusters_indices = []
        for line in gene_positions:
            gene_index = line
            for row in condition_positions:
                each_bicluster = []
                each_bicluster_index = []
                for item in gene_index:
                    each_row = []
                    for column_index in row:
                        each_row.append(gene_normalized_array[item,column_index])
                    each_bicluster.append(each_row)
                    each_bicluster_index.append(gene_index)
                    each_bicluster_index.append(row)
                biclusters.append(each_bicluster)
                biclusters_indices.append(each_bicluster_index)
        return biclusters, biclusters_indices

    #MSR function for each bicluster
    def MSR_function(self, biclusters):
        biclusters_MSR = []
        for line in biclusters:
            bicluster_residue = []
            aiJ_list = []
            aIj_list = []

            #calculate the mean of each row in the bicluster
            for row in line:
                row_mean = mean(row)
                aiJ_list.append(row_mean)

            #calculate the mean of each column in the bicluster
            bicluster_transpose = np.array(line).T.tolist()
            for column in bicluster_transpose:
                column_mean = mean(column)
                aIj_list.append(column_mean)

            aIJ = mean(aiJ_list)    #mean of the bicluster

            #calculate the residue of each element in the bicluster
            for i, row in enumerate(line):
                row_residue = []
                for j, item in enumerate(row):
                    res = item + aIJ - aiJ_list[i] - aIj_list[j]
                    row_residue.append(res)
                bicluster_residue.append(row_residue)

            #calculate the MSR of the bicuster
            residue_sum = 0.0
            I = len(bicluster_residue)
            J = len(bicluster_residue[0])

            for row in bicluster_residue:
                row_sum = 0.0
                for item in row:
                    row_sum = row_sum + item**2
                residue_sum = residue_sum + row_sum
            MSR = residue_sum/(I*J)
            biclusters_MSR.append(MSR)
        return biclusters_MSR

    #function that removes biclusters with MSR>threshold(yeast=300, schisto=2100, lymphoma=1200)
    def bicluster_seeding(self, biclusters_MSR, biclusters, biclusters_indices):
        for j, item in enumerate(biclusters_MSR):
            if biclusters_MSR[j]>2100:
                del(biclusters[j])
                del(biclusters_indices[j])
                del(biclusters_MSR[j])
        return biclusters_MSR, biclusters, biclusters_indices

    # function for ENCODING
    def encoding(self, biclusters_indices):
        data_dimen = gene_normalized_array.shape
        encoded_biclusters = []
        result_file = open("encoded.csv",'w')
        wr = csv.writer(result_file, dialect='excel')
        for line in biclusters_indices:
            g_index = line[0]
            c_index = line[1]
            gene_code = np.zeros(data_dimen[0], dtype = int)
            condition_code = np.zeros(data_dimen[1], dtype = int)
            for item in g_index:
                gene_code[item] = 1  
            for item in c_index:
                condition_code[item] = 1
            encoded_bi = gene_code.tolist() + condition_code.tolist()
            encoded_biclusters.append(encoded_bi)
            wr.writerow(encoded_bi)
        return encoded_biclusters

