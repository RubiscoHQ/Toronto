from random import sample
import pandas as pd
from pandas import DataFrame
from numpy import array
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn import linear_model
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import time
import os
import json


# File routes
import_dir = '/Users/rubisco/Desktop/Toronto/Mutation_analyse/input/'
output_dir = '/Users/rubisco/Desktop/Toronto/Mutation_analyse/output_%s/'\
             % time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime())
patho_file = 'clinvar.hg19_multianno.filt.new.narm.txt'
ctrl_file = 'ExAC.all.hg19_multianno.new.vcf.hq.input.narm.txt'
out_file_name = 'analysis.log.txt'

# Score to analyze
score_list = ['SIFT']
group_select_list = ['I', 'II', 'III', 'IV', 'V', 'VI']

# Contaminate setting
contaminate_select = 9000
contaminate_start = 0.00
contaminate_end = 0.40
contaminate_step = 0.01

os.mkdir(output_dir)
out_file = open(output_dir + out_file_name, 'w')
column_name = ['Chr', 'Pos', 'Ref', 'Alt', 'Gene_refGene',
               'SIFT', 'Polyphen2_HDIV', 'Polyphen2_HVAR', 'MutationTaster', 'MutationAssessor',
               'PROVEAN', 'VEST3', 'CADD', 'DANN', 'fathmm-MKL_coding', 'MetaSVM', 'MetaLR',
               'GERP++_RS', 'phyloP7way_vertebrate', 'phyloP20way_mammalian',
               'phastCons7way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds',
               'AC', 'AN', 'AF', 'Freq_group']
out_format = ['prior', 'pathogenic_mutation', 'control_mutation', 'score', 'freq_group',
              'select', 'contaminate', 'contaminate_rate', 'auprc', 'pr90', 'pr10']
f = open(output_dir + 'settings.log.txt', 'w')
setting = {'patho': patho_file, 'ctrl': ctrl_file, 'score': score_list, 'group': group_select_list,
           'contaminate': {'select': contaminate_select, 'start': contaminate_start,
                           'end': contaminate_end, 'step': contaminate_step}}
f.write(json.dumps(setting, indent=4))


def split_df_gene(df):
    td = {}
    value = df.values
    for line in value:
        if line[4] not in td:
            td[line[4]] = [list(line)]
        else:
            td[line[4]].append(list(line))
    return td


def get_score_df(patho, ctrl, score, gene_spec=False, random=True, freqg=None, contaminate=None):
    if freqg:
        ctrl = ctrl[ctrl['Freq_group'] == freqg]
        if len(ctrl) == 0:
            raise MyError('No freq group')

    if gene_spec:
        total_plist = []
        total_clist = []
        patho_df_dict = split_df_gene(patho)
        ctrl_df_dict = split_df_gene(ctrl)
        for gene, pvalue in patho_df_dict.items():
            if gene not in ctrl_df_dict:
                continue
            patho_count = len(pvalue)
            cvalue = ctrl_df_dict[gene]
            ctrl_count = len(cvalue)
            if patho_count > ctrl_count:
                pvalue = sample(pvalue, ctrl_count)
            else:
                cvalue = sample(cvalue, len(pvalue))
            total_plist += pvalue
            total_clist += cvalue

        patho = DataFrame(total_plist, columns=column_name)
        ctrl = DataFrame(total_clist, columns=column_name)
    else:
        if not contaminate:
            if random:
                if len(patho) > len(ctrl):
                    patho = patho.sample(n=len(ctrl))
                else:
                    ctrl = ctrl.sample(n=len(patho))
        else:
            try:
                patho = patho.sample(n=int(contaminate_select + contaminate_select * contaminate))
            except ValueError:
                print 'Pathogenic varitants not enough.', 'Contaminate rate is', contaminate, \
                '\nNeeds', int(contaminate_select + contaminate_select * contaminate), 'Have', len(patho)
                print 'Not perform contamination.'
                raise

            pieces = [patho[:contaminate_select], patho[contaminate_select:]]
            patho = pieces[0]

            try:
                ctrl = pd.concat([ctrl.sample(n=int(contaminate_select - contaminate_select * contaminate)),
                                  pieces[1]])
            except ValueError:
                print 'ExAC varitants not enough.', 'Contaminate rate is', contaminate, \
                '\nNeeds', int(contaminate_select - contaminate_select * contaminate), 'Have', len(ctrl)
                print 'Not perform contamination.'
                raise

    # Information selection
    score_df = DataFrame({'score': list(patho[score]) + list(ctrl[score]),
                          'logic': [True] * len(patho) + [False] * len(ctrl),
                          'Gene_refGene': list(patho['Gene_refGene']) + list(ctrl['Gene_refGene']),
                          'AF': list(patho['AF']) + list(ctrl['AF']),
                          'Freq_group': list(patho['Freq_group']) + list(ctrl['Freq_group'])
                          })
    out_file.write('%.3f' % (float(len(patho)) / (len(ctrl) + len(patho))) +
                     '\t' + '%d' % len(patho) + '\t' + '%d' % len(ctrl) + '\t')
    if score_df is None:
        raise
    return score_df


def calc_PRexat(precision, recall, exat):
    precision = list(precision)
    recall = list(recall)
    exat = float(exat)
    if not (0 <= exat <= 1):
        raise
    upexat = 1-exat
    up = 1.0
    downexat = -exat
    down = 0.0
    for i in recall:
        if 0 < i-exat < upexat:
            up = i
            upexat = i-exat
        elif downexat < i-exat < 0:
            down = i
            downexat = i-exat
    up_value = upexat/(upexat-downexat)
    down_value = 1-up_value
    up_precision = precision[recall.index(up)]
    down_precision = precision[recall.index(down)]
    prexat = up_precision * down_value + down_precision * up_value  # the closer the more value, lever principle
    return prexat


class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def neg(value):
    return 1-value


def process_score_df(score_df):
    precision, recall, thresholds = precision_recall_curve(array(score_df['logic']), array(score_df['score']))
    auprc = average_precision_score(array(score_df['logic']), array(score_df['score']))
    if auprc < 0.5:
        precision = map(neg, precision)
        auprc = neg(auprc)
    pr90 = calc_PRexat(precision, recall, 0.90)
    pr10 = calc_PRexat(precision, recall, 0.10)
    pgene_count = score_df[score_df['logic'] == True].groupby(['Gene_refGene'])['AF'].count().to_dict()
    cgene_count = score_df[score_df['logic'] == False].groupby(['Gene_refGene'])['AF'].count().to_dict()
    pmut_freq = score_df[score_df['logic'] == True]['AF']
    cmut_freq = score_df[score_df['logic'] == False]['AF']
    pscore = score_df[score_df['logic'] == True]['score']
    cscore = score_df[score_df['logic'] == False]['score']
    return precision, recall, auprc, pr90, pr10, pgene_count, cgene_count, pmut_freq, cmut_freq, pscore, cscore


def analysis_data():
    # Read input
    patho_df = pd.read_table(import_dir + patho_file)
    ctrl_df = pd.read_table(import_dir + ctrl_file)

    # Setup output
    out_file.write('\t'.join(out_format)+'\n')
    td = {}

    # Calculate
    for score in score_list:
        print score
        for group in group_select_list:
            print group,
            try:
                for select in ['genome', 'gss']:
                    print select,
                    if select == 'genome':  # Genome wide select control variants
                        for contaminate in ['T', 'F']:
                            print contaminate,
                            if contaminate == 'F':
                                # No contamination
                                index = (score, group, select, 'F', '0')
                                score_df = get_score_df(patho=patho_df, ctrl=ctrl_df, gene_spec=False,
                                                        score=score, random=True, freqg=group, contaminate=None)
                                out_file.write('\t'.join(index))
                                td[index] = process_score_df(score_df)
                                out_file.write('\t' + '\t'.join(map(str, td[index][2:5])) + '\n')

                            else:
                                # contamination
                                c = contaminate_start
                                while c <= contaminate_end:
                                    index = (score, group, select, contaminate, '%.2f' % c)
                                    score_df = get_score_df(patho=patho_df, ctrl=ctrl_df, gene_spec=False,
                                                            score=score, random=True, freqg=group, contaminate=c)
                                    out_file.write('\t'.join(index))
                                    td[index] = process_score_df(score_df)
                                    out_file.write('\t' + '\t'.join(map(str, td[index][2:5])) + '\n')
                                    c += contaminate_step

                    else:  # gene specific select control variants
                        index = (score, group, select, 'F', '0')
                        score_df = get_score_df(patho=patho_df, ctrl=ctrl_df, gene_spec=True,
                                                score=score, random=True, freqg=group, contaminate=None)
                        out_file.write('\t'.join(index))
                        td[index] = process_score_df(score_df)
                        out_file.write('\t' + '\t'.join(map(str, td[index][2:5])) + '\n')
            except MyError:
                continue
            print
        print
    return td


def draw_prc(data_dict, prefix, contaminate=False):
    plt.clf()
    items = data_dict.items()
    items.sort()
    for key, value in items:
        plt.plot(value[1], value[0], label='{2} AUC={0:0.3f} PR90={1:0.3f}'
                 .format(value[2], value[3], key))

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title(prefix + ' Precision-Recall Curve')
    if not contaminate:
        plt.legend(loc="lower left")
    plt.savefig(output_dir + prefix + '.prc.png')
    return


def draw_score_distribution(data_dict, prefix, sel=1, shade=True):
    plt.clf()
    items = data_dict.items()
    items.sort()
    for key, value in items:
        plot = sns.distplot(array(value[sel]), hist=False, rug=False, kde_kws={"shade": shade}, label=key)
        plt.xlabel('Score value')
        plt.title(prefix + ' Score Distribution')
        plt.legend(loc="upper left")
    named = {1: 'control', 0: 'patho'}
    plot.get_figure().savefig(output_dir + prefix + '_' + named[sel] +'.dist.png')
    return


def draw_gene_count(data_dict, prefix):
    for key, value in data_dict.items():
        plt.clf()
        pd = value[0]
        cd = value[1]
        td = {}
        for gene in pd:
            if gene in cd:
                td[gene] = [pd[gene], cd[gene]]
            else:
                td[gene] = [pd[gene], 0]
        for gene in cd:
            if gene in td:
                continue
            else:
                td[gene] = [0, cd[gene]]
        x = []
        y = []
        for gkey, gvalue in td.items():
            x.append(gvalue[0])
            y.append(gvalue[1])
        sns_plot = sns.jointplot(array(x), array(y), stat_func=None) \
            .plot_marginals(sns.distplot, kde=False, color=".5")
        sns_plot.savefig(output_dir + prefix + '_' + key + '.count.png')


def draw_plot(td):
    for score in score_list:
        #  Draw plot for frequent group prepare
        freq_pr_dict = {}
        freq_score_dict = {}
        freq_gene_dict = {}
        prefix = '_'.join([score, 'Frequent'])
        for group in ['I', 'II', 'III', 'IV', 'V', 'VI']:
            for key, value in td.items():
                if (score, group, 'genome', 'F') != key[:4]:
                    continue
                freq_pr_dict[key[1]] = value[:5]
                freq_score_dict[key[1]] = value[-2:]
                freq_gene_dict[key[1]] = value[5:7]
        draw_prc(freq_pr_dict, prefix)  # freq group PRC
        draw_score_distribution(freq_score_dict, prefix, sel=1)  # ctrl mut dist against selection
        draw_score_distribution(freq_score_dict, prefix, sel=0)  # patho mut dist against selection
        #  draw_gene_count(freq_gene_dict, prefix)  # gene count plot in ctrl and patho

        #  Draw plot for contaminant rate
        freq_pr_contaminate_dict = {}
        freq_score_contaminate_dict = {}
        for key, value in td.items():
            if score != key[0] or ('genome', 'T') != key[2:4]:
                continue
            group = key[1]
            if group not in freq_score_contaminate_dict:
                freq_score_contaminate_dict[group] = {key[4]: value[-2:]}
                freq_pr_contaminate_dict[group] = {key[4]: value[:5]}
            else:
                freq_score_contaminate_dict[group][key[4]] = value[-2:]
                freq_pr_contaminate_dict[group][key[4]] = value[:5]

        for group, gdict in freq_score_contaminate_dict.items():
            prefix = '_'.join([score, group, 'contamination'])
            #  draw_score_distribution(gdict, prefix, shade=False)  # ctrl mut dist against

        for group, gdict in freq_pr_contaminate_dict.items():
            prefix = '_'.join([score, group, 'contamination'])
            draw_prc(gdict, prefix, contaminate=True)  # contaminate rate PRC

    return


def analysis_prc():
    data0 = pd.read_table(output_dir + out_file_name)

    # For contamination analysis
    data = data0[data0.contaminate == 'T']
    for score in list(set(list(data.score))):
        data = data[data.score == score]

        func_dict = {}
        data_dict = {}
        group_list = list(set(list(data.freq_group)))
        group_list.sort()
        for group in group_list:
            data_group = data[data.freq_group == group]
            fit = linear_model.LinearRegression()
            b = array(data_group.auprc).transpose()
            a = array([data_group.contaminate_rate]).transpose()
            y = fit.fit(a, b).predict(a)
            func_dict[group] = fit
            data_dict[group] = [a, b, y]

        grey = 0.1
        plt.clf()
        for group, data_l in data_dict.items():
            plt.scatter(data_l[0], data_l[1], color=str(grey))
            plt.plot(data_l[0], data_l[2], color=str(grey), linewidth=2)
            grey += 0.15
        plt.xlim([-0.01, contaminate_end])
        plt.ylim([min(data.auprc) - 0.003, max(data.auprc) + 0.003])
        plt.xlabel('Contamination Rate')
        plt.ylabel('AUPRC')
        plt.legend(group_list)
        plt.savefig(output_dir + score + '_contaminate_scatter.png')

        x = 0.0
        predict_dict = {}
        predict_dict_0 = {}
        while x < 1.0:
            for group, func in func_dict.items():
                if group not in predict_dict:
                    predict_dict_0[group] = func.predict(0.0)[0]
                    predict_dict[group] = [[func.predict(x)[0], x]]
                else:
                    predict_dict[group].append([func.predict(x)[0], x])
            x += 0.001

        result_dict = {}
        items = predict_dict.items()
        items.sort()
        for group, values in items:
            result_dict[group] = []
            items0 = predict_dict_0.items()
            items0.sort()
            for group_0, value_0 in items0:
                larger = []
                larger_x = []
                smaller = []
                smaller_x = []
                for value in values:
                    if value[0] < value_0:
                        smaller.append(value[0])
                        smaller_x.append(value[1])

                    else:
                        larger.append(value[0])
                        larger_x.append(value[1])
                if len(larger) == 0:
                    result_dict[group].append(0.0)
                    continue
                larger_value = min(larger)
                smaller_value = max(smaller)
                larger_value_x = larger_x[larger.index(larger_value)]
                smaller_value_x = smaller_x[smaller.index(smaller_value)]
                larger_large = larger_value-value_0
                smaller_small = value_0 - smaller_value
                if larger_large <= smaller_small:
                    result_dict[group].append(larger_value_x)
                else:
                    result_dict[group].append(smaller_value_x)

        rd = pd.DataFrame(result_dict, index=group_list)
        plt.clf()
        mask = np.zeros_like(rd.T)
        mask[np.triu_indices_from(mask)] = True
        with sns.axes_style("white"):
            sns.heatmap(rd.T, linewidths=1, annot=True, mask=mask)
            plt.title(score+' Relative Contamination Rate')
            plt.savefig(output_dir + score+'_Relative_Contamination_Rate.count.png')

    return

if __name__ == '__main__':
    td = analysis_data()
    draw_plot(td)
    out_file.close()
    analysis_prc()

