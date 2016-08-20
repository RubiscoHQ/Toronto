titlel = ['Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene',
          'AAChange.refGene','CLINSIG','CLNDBN','CLNACC','CLNDSDB','CLNDSDBID',
          'ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','ExAC_nonpsych_ALL',
          'ExAC_nonpsych_AFR','ExAC_nonpsych_AMR','ExAC_nonpsych_EAS','ExAC_nonpsych_FIN','ExAC_nonpsych_NFE','ExAC_nonpsych_OTH','ExAC_nonpsych_SAS','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_pred','LRT_score','LRT_pred','MutationTaster_score','MutationTaster_pred','MutationAssessor_score','MutationAssessor_pred','FATHMM_score','FATHMM_pred','PROVEAN_score','PROVEAN_pred','VEST3_score','CADD_raw','CADD_phred','DANN_score','fathmm-MKL_coding_score','fathmm-MKL_coding_pred','MetaSVM_score','MetaSVM_pred','MetaLR_score','MetaLR_pred','integrated_fitCons_score','integrated_confidence_value','GERP++_RS','phyloP7way_vertebrate','phyloP20way_mammalian','phastCons7way_vertebrate','phastCons20way_mammalian','SiPhy_29way_logOdds', 'QUAL','FILTER','INFO']

psy = '/Volumes/lhq_data/data/db/ExAC.nonpsy.hg19_multianno.txt'
tcga = '/Volumes/lhq_data/data/db/ExAC.nontcga.hg19_multianno.txt'
nf = '/Volumes/lhq_data/data/db/ExAC.nonpat.hg19_multianno.txt'


def select(anno):
    f = open(anno)
    nfanno = open(anno[:-4] + '.new.txt', 'w')
    nfannop = open(anno[:-4] + '.patho.txt', 'w')

    title = True

    for line in f:
        if title:
            title = False
            nfanno.write(line)
            continue
        ll = line.strip().split('\t')

        # 'ExonicFunc.refGene'
        if ll[titlel.index('ExonicFunc.refGene')] != 'nonsynonymous SNV':
            continue

        if 'Pathogenic' in ll[-5].split('|'):
            nfannop.write(line)
            continue


        nfanno.write(line)
    return

# select('/Volumes/lhq_data/data/db/ExAC.nontcga.hg19_multianno.txt')
# select('/Volumes/lhq_data/data/db/ExAC.nonpst.hg19_multianno.txt')
# select('/Volumes/lhq_data/data/db/ExAC.nonpat.hg19_multianno.txt')
# select('/Volumes/lhq_data/data/db/ExAC.all.hg19_multianno.txt')



def merge_var(r1, r2, r3):
    f1 = open(r1)
    td = {}
    for line in f1:
        key = tuple(line.strip().split()[:5])
        if key in td:
            print key
        else:
            td[key] = line
    f1.close()
    print len(td)

    f2 = open(r2)
    nf = open(r3, 'w')
    n=0
    for line in f2:
        key = tuple(line.strip().split()[:5])
        if key in td:
            nf.write(td[key])
            del td[key]
            n+=1
            if n%100000==0:
                print n,

def anno_exac_vcf(anno, vcf):
    fanno = open(anno)
    nfanno = open(anno[:-4]+'.vcf.txt', 'w')
    ufanno = open(anno[:-4] + '.novcf.txt', 'w')
    title = True
    td = {}
    for line in fanno:
        if title:
            title = False
            nfanno.write(line)
            continue
        ll = line.strip().split('\t')
        key = tuple(ll[:2] + ll[3:5])
        if key in td:
            print key
            raise
        td[key] = ll

    print len(td), 'total annoed var'

    n = 0
    fvcf = open(vcf)
    for line in fvcf:
        if line[0] == '#':
            continue
        ll = line.strip().split('\t')
        key = tuple(ll[:2] + ll[3:5])
        if key in td:
            nfanno.write('\t'.join(td[key]+ll[5:])+'\n')
            del td[key]
            n+=1
            if n%100000 == 0:
                print n,

    print len(td), 'unvcf var'
    for key in td:
        ufanno.write('\t'.join(td[key]) + '\n')
# anno_exac_vcf('/Volumes/lhq_data/data/db/ExAC.nonpat.hg19_multianno.new.txt', '/Volumes/lhq_data/data/db/ExAC.r0.3.1.sites.vep.vcf')

def exac_hq(r, onehom=False):
    f = open(r)
    if onehom:
        nf = open(r[:-4]+'.hq.1hom.txt', 'w')
        no_nf = open(r[:-4]+'.hq.not1hom.txt', 'w')
    else:
        nf = open(r[:-4]+'.hq.txt', 'w')
    freq_region = {14: [[106329000, 106331000], [107178000, 107180000]],
                   2: [[89160000, 89162000]],
                   17: [[18967000, 18968000], [19091000, 19092000]],
                   22: [[23223000, 23224000]],
                   1: [[152975000, 152976000]],
                   }
    title = True
    n = 0
    for line in f:
        if title:
            title = False
            nf.write(line)
            continue
        n += 1
        if n % 100000 == 0:
            print n
        ll = line.strip().split('\t')
        if ll[-2] != 'PASS':
            continue
        infol = ll[-1].split(';')
        infod = {}
        for eq in infol:
            eql = eq.split('=')
            try:
                infod[eql[0]] = eql[1]
            except IndexError:
                infod[eql[0]] = 'NA'

        if int(infod['AC_Adj']) < 1:
            continue
        if int(infod['AN_Adj']) < 97130:
            continue

        try:
            chrom = int(ll[0])
            if chrom in freq_region:
                for pair in freq_region[chrom]:
                    if pair[0] < int(ll[1]) < pair[1]:
                        continue
        except ValueError:
            pass

        if onehom:
            if int(infod['AC_Hom']) < 1:
                no_nf.write(line)
                continue

        nf.write(line)

#exac_hq('/Volumes/lhq_data/data/db/ExAC.nonpat.hg19_multianno.new.vcf.txt', True)


def retrench_file(r, s):
    f = open(r)
    nf = open(s, 'w')
    title = True
    for line in f:
        if title:
            title = False
            tl = ['Chr','Pos','Ref','Alt','Gene_refGene','SIFT','Polyphen2_HDIV',
                  'Polyphen2_HVAR','LRT','MutationTaster','MutationAssessor','FATHMM',
                  'PROVEAN','VEST3','CADD','DANN','fathmm-MKL_coding','MetaSVM',
                  'MetaLR','GERP++_RS','phyloP7way_vertebrate','phyloP20way_mammalian',
                  'phastCons7way_vertebrate','phastCons20way_mammalian','SiPhy_29way_logOdds',
                  'AC','AN','AF','Freq_group']
            nf.write('\t'.join(tl)+'\n')
            continue
        nll = []
        ll = line.strip().split('\t')
        nll = nll + ll[:1] + ll[2:5]
        nll.append(ll[6])  # gene

        scores = [ll[10], ll[12], ll[14], ll[16], ll[18], ll[20], ll[22], ll[24], ll[26], ll[28],
                  ll[29], ll[30], ll[32], ll[34], ll[38], ll[39], ll[40], ll[41], ll[42], ll[43]]
        new_scores = []
        for score in scores:
            if score in ['.', 'NA']:
                new_scores.append('NA')
            else:
                new_scores.append(score)


        nll += new_scores  # scores
        infol = ll[-1].split(';')
        infod = {}
        for eq in infol:
            eql = eq.split('=')
            try:
                infod[eql[0]] = eql[1]
            except IndexError:
                infod[eql[0]] = 'NA'
        nll += [infod['AC'], infod['AN'], infod['AF']]

        #  Freq group
        if int(infod['AC']) == 1:
            nll.append('I')
        elif 2 <= int(infod['AC']) <= 12:
            nll.append('II')
        elif 0.0001 <= float(infod['AF']) <= 0.001:
            nll.append('III')
        elif 0.001 <= float(infod['AF']) <= 0.01:
            nll.append('IV')
        elif 0.01 <= float(infod['AF']) <= 0.1:
            nll.append('V')
        else:
            nll.append('VI')

        nf.write('\t'.join(nll)+'\n')
    return


#retrench_file('/Volumes/lhq_data/data/db/ExAC.all.hg19_multianno.new.vcf.hq.1hom.txt', '/Volumes/lhq_data/data/db/ExAC.all.hg19_multianno.new.vcf.hq.1hom.input.txt')


def change_clinvar(r):
    f = open(r)
    nf = open(r[:-4]+'.new.txt', 'w')
    title = True
    for line in f:
        if title:
            title = False
            tl = ['Chr', 'Pos', 'Ref', 'Alt', 'Gene_refGene', 'SIFT', 'Polyphen2_HDIV',
                  'Polyphen2_HVAR', 'LRT', 'MutationTaster', 'MutationAssessor', 'FATHMM',
                  'PROVEAN', 'VEST3', 'CADD', 'DANN', 'fathmm-MKL_coding', 'MetaSVM',
                  'MetaLR', 'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian',
                  'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds',
                  'AC', 'AN', 'AF', 'Freq_group']
            nf.write('\t'.join(tl) + '\n')
            continue
        nll = []
        ll = line.strip().split('\t')
        nll = nll + ll[:1] + ll[2:5]
        nll.append(ll[6])  # gene
        scores = [ll[31], ll[33], ll[35], ll[37], ll[39], ll[41], ll[43], ll[45], ll[47], ll[49],
                  ll[50], ll[51], ll[53], ll[55], ll[59], ll[60], ll[61], ll[62], ll[63], ll[64]]
        new_scores = []
        for score in scores:
            if score in ['.', 'NA']:
                new_scores.append('NA')
            else:
                new_scores.append(score)
        nll += new_scores  # scores
        nll += ['NA', 'NA', ll[15]]
        nll.append('P')  # Patho group
        nf.write('\t'.join(nll) + '\n')
    return


def input_narm(r):
    f = open(r)
    nf = open(r[:-4]+'.narm.txt', 'w')

    for line in f:
        ll = line.strip().split('\t')
        score = ll[5:8]+ll[9:11]+ll[12:25]
        flag = False
        for i in score:
            if i == 'NA':
                flag = True
        if flag:
            continue
        nll = ll[:5]+score+ll[-4:]
        nf.write('\t'.join(nll)+'\n')

input_narm('/Volumes/lhq_data/data/db/ExAC.nonpat.hg19_multianno.new.vcf.hq.input.txt')


def select_group(file, group):
    f = open(file)
    nf = open(file[:-4]+'.'+group+'.txt', 'w')
    for line in f:
        ll = line.strip().split('\t')
        if ll[-1] == group:
            nf.write(line)

#select_group('/Volumes/lhq_data/data/db/ExAC.nonpat.hg19_multianno.new.vcf.hq.not1hom.input.txt', 'VI')