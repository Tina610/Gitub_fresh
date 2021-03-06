#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/6/1 下午3:58
# @Author  : liting
# @contact : tinating610@163.com
# @File    : chemmodula_V1.py
# @Software: PyCharm

import os, sys
import re
import Dict
import xlrd
import itertools
from argparse import ArgumentParser


class CHEMdrug():
    def __init__(self, chemhap, drug, multimutaion, c1info):
        self.haplotype = chemhap
        self.drug = drug
        self.c1info = c1info
        self.annofile = multimutaion
        self.haplotypedict = {}
        self.drugdict = {}
        self.endchesite = {}
        global rslist
    
    def readhaptype(self):
        global rslist
        rslist = []
        work = xlrd.open_workbook(self.haplotype)
        sheet = work.sheet_by_index(0)
        nrows = sheet.nrows
        for row in range(1, nrows):
            rowvalue = sheet.row_values(row)
            drug = rowvalue[0]
            gene = rowvalue[1]
            genetype = rowvalue[2]
            clinal = str(rowvalue[5]).replace('%', r'\%')
            conclusion = str(rowvalue[6]).replace('%', r'\%')
            if sheet.cell(row, 3).ctype == 0:
                basic = '-'
            elif sheet.cell(row, 3).value == 'others':
                basic = '无'
            else:
                basic = rowvalue[3]
            Dict.addtodict3(self.haplotypedict, drug, gene, genetype, [basic, clinal, conclusion])
            rslist.append(rowvalue[3])
    
    def readdrug(self):
        work = xlrd.open_workbook(self.drug)
        sheet = work.sheet_by_index(0)
        nrows = sheet.nrows
        for row in range(1, nrows):
            rowvalue = sheet.row_values(row)
            if sheet.cell(row, 4).ctype == 0:
                ishaplotype = 'N'
            else:
                ishaplotype = 'H'
            Dict.addtodict5(self.drugdict, rowvalue[0], rowvalue[1], rowvalue[2], rowvalue[3], ishaplotype, rowvalue[5])
    
    def readannofile(self):
        global rslist
        mutationrs = {}
        with open(self.annofile, 'r') as F:
            mark = None
            title = F.readline().strip().split('\t')
            for idx, value in enumerate(title):
                if 'snp13' in value:
                    mark = idx
            
            for line in F:
                lines = line.strip('\n').split('\t')
                if lines[mark] in rslist:
                    if lines[mark] == 'rs2032582' and re.search(r'^1/2:', lines[-1]):
                        Dict.addtodict(mutationrs, lines[mark], 'CT')
                    else:
                        gt = self.percent(lines)
                        Dict.addtodict(mutationrs, lines[mark], gt)
        
        return mutationrs
    
    def searchdrug(self):
        mutationrs = self.readannofile()
        for drug in self.drugdict.keys():
            for gene in self.drugdict[drug].keys():
                rscom = []
                rssite = {}
                for cloc in self.drugdict[drug][gene].keys():
                    for rsid in self.drugdict[drug][gene][cloc].keys():
                        for hap in self.drugdict[drug][gene][cloc][rsid].keys():
                            genechem = gene + '(' + cloc + ')'
                            genetype = ''
                            if hap == 'N':
                                if rsid in mutationrs:
                                    genetype = rsid + ':' + mutationrs[rsid]
                                else:
                                    genetype = rsid + ':' + self.drugdict[drug][gene][cloc][rsid][hap]
                                
                                if genetype in self.haplotypedict[drug][gene].keys():
                                    Dict.addtodict3(self.endchesite, drug, 'N', genechem,
                                                    [genetype.split(':')[1]] + self.haplotypedict[drug][gene][genetype])
                                else:
                                    print('strange genetype: ' + hap + ' ' + rsid + ' ' + gene + drug)
                            else:
                                if rsid in mutationrs:
                                    genetype = rsid + ':' + mutationrs[rsid]
                                    rssite[genetype] = [genechem, mutationrs[rsid]]
                                else:
                                    genetype = rsid + ':' + self.drugdict[drug][gene][cloc][rsid][hap]
                                    rssite[genetype] = [genechem, self.drugdict[drug][gene][cloc][rsid][hap]]
                                rscom.append(genetype)
                if rscom:
                    combiners = self.comiter(rscom)
                    mark = 0
                    for coms in combiners:
                        if coms in self.haplotypedict[drug][gene].keys():
                            mark += 1
                            for site in coms.split(','):
                                if site in rssite:
                                    Dict.addtodict3(self.endchesite, drug, 'H', rssite[site][0], [rssite[site][1]] +
                                                    self.haplotypedict[drug][gene][coms])
                                else:
                                    print('this rssite is not in list ' + site)
                        else:
                            pass
                    if mark == 0:
                        for site in rssite:
                            Dict.addtodict3(self.endchesite, drug, 'H', rssite[site][0],
                                            [rssite[site][1]] + self.haplotypedict[drug][gene]['others'])
                    else:
                        print(gene + ' is haplotype')
    
    def comiter(self, listscom):
        combine = []
        for i in range(2, 7):
            for com in itertools.permutations(listscom, i):
                tmp = ','.join(com)
                combine.append(tmp)
        return combine
    
    def zh_chinese(self, strs):
        zh_pattern = re.compile(u'[\u4e00-\u9fff]+')
        if re.findall(zh_pattern, strs):
            return True
        else:
            return False
    
    def numlin(self, strs, one):
        new = ''
        j = 0
        lines = []
        for i in strs:
            if self.zh_chinese(i):
                j += 1  ###3.36
                new = new + i
            elif i.isdigit():
                new = new + i
                j += 0.48
            elif i.islower():
                new = new + i
                j += 0.5
            elif i.isupper():
                new = new + i
                j += 0.65  ###2.36
            elif re.match(r'，|。|%|、|；', i):
                new = new + i
                j += 0.8
            else:
                new = new + i
            if j > one:
                # new = new + r'\newline '
                lines.append(new)
                new = ''
                j = 0
            else:
                pass
        nums = len(lines) + 1
        return nums
    
    def computlines(self):
        temp_chemsit = {}
        all_num_tex = {}
        for key1 in self.endchesite.keys():
            drug_lines = 0
            for key2 in self.endchesite[key1].keys():
                if key2 == 'N':
                    for key3 in self.endchesite[key1][key2].keys():
                        [basic, hploy, clinal, conclusion] = self.endchesite[key1][key2][key3]
                        n_clinal = self.numlin(clinal, 35)
                        n_conclusion = self.numlin(conclusion, 3)
                        Dict.addtodict3(temp_chemsit, key1, key3, 'conclusion', n_conclusion)
                        Dict.addtodict3(temp_chemsit, key1, key3, 'clinal', n_clinal)
                        drug_lines = drug_lines + n_clinal * 0.7 + n_conclusion
                else:
                    key3list = list(self.endchesite[key1][key2].keys())
                    for key3 in self.endchesite[key1][key2].keys():
                        [basic, hploy, clinal, conclusion] = self.endchesite[key1][key2][key3list[0]]
                        n_clinal = self.numlin(clinal, 35)
                        Dict.addtodict3(temp_chemsit, key1, key3, 'conclusion', len(key3list))
                        Dict.addtodict3(temp_chemsit, key1, key3, 'clinal', n_clinal)
                    drug_lines = drug_lines + len(key3list) + n_clinal * 0.7
            all_num_tex[key1] = int(drug_lines)
        return temp_chemsit, all_num_tex
    
    def outcheck(self, checkout):
        with open(checkout, 'w') as OUT:
            for drug in self.endchesite:
                for hap in self.endchesite[drug].keys():
                    for site in self.endchesite[drug][hap].keys():
                        gene, cloc = re.search('(\w+)\((..+)\)', site).groups()
                        OUT.write(gene + '\t' + cloc + '\t' + self.endchesite[drug][hap][site][0] + '\n')
    
    def modifyC1(self, out):
        temp_chemsite, all_num_tex = self.computlines()
        OUT = open(out, 'w')
        total = 7
        init_len = 0
        with open(self.c1info, 'r') as F:
            for line in F:
                print(line, init_len, total)
                for drug in sorted(self.endchesite.keys()):
                    for hap in sorted(self.endchesite[drug].keys()):
                        for site in sorted(self.endchesite[drug][hap].keys()):
                            gene = str(site).split('(')[0]
                            allnums = all_num_tex[drug]
                            # print(init_len,total,'this is len init_len')
                            if re.search('{numline}{\*}' + '{' + drug + '}', line):
                                print(init_len, 'numline++++++', drug, allnums)
                                if init_len > total:
                                    # print('outline +++++++++++++++++++++++++----------')
                                    init_len = 0
                                    total = 34
                                    # line=re.sub(r'cline{\d+-\d+}','hline',line)
                                    OUT.write('\end{xltabular}\n')
                                    OUT.write(r'\newpage' + '\n')
                                    OUT.write(
                                        r'\begin{xltabular}{\textwidth}{| m{1.8cm}<{\centering} | m{3.5cm}<{\centering} | m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.7cm}<{\centering}| m{2cm}<{\centering} | m{0.8cm}<{\centering} |}' + '\n'
                                        + r'\arrayrulecolor{深蓝}\hline 药物 & 检测位点 & 检测结果 & 单体型 & 多态性类型 & 结论 & 等级 \\\hline' + '\n')
                                if drug == '6-巯基嘌呤': allnums = 4
                                line = re.sub('numline', str(allnums), line)
                                init_len += temp_chemsite[drug][site]['conclusion']
                                print(init_len, drug, site)
                            if re.search(r'label.haptype.' + drug + '.' + gene, line):
                                # print('label.haptype.',line,gene)
                                line = re.sub(r'label.haptype.' + drug + '.' + gene,
                                              self.endchesite[drug]['H'][site][1], line)
                            if re.search(r'label.toxicity.' + drug + '.' + gene, line):
                                # print('label.toxicity.')
                                line = re.sub(r'label.toxicity.' + drug + '.' + gene,
                                              self.endchesite[drug]['H'][site][3], line)
                            # if re.search(site, line) and re.search(r'XX', line):
                            if line.find(site) != -1 and line.find('XX') != -1:
                                init_len += 1
                                if init_len > total:
                                    print(init_len, total, line)
                                    # init_len = 0
                                    # total=34
                                    line = re.sub(r'cline{2-7}', r'hline', line)
                                line = re.sub(r'XX', self.endchesite[drug][hap][site][0], line)
                                # init_len += 1
                            if re.search(r'label.haptype.' + drug + '.clin.' + gene, line):
                                if init_len > total:
                                    print(init_len, total, 'label.haptype.=====cline', line)
                                    init_len = 0
                                    total = 34
                                    OUT.write('\end{xltabular}\n')
                                    OUT.write(r'\newpage' + '\n')
                                    OUT.write(
                                        r'\begin{xltabular}{\textwidth}{| m{1.8cm}<{\centering} | m{3.5cm}<{\centering} | m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.7cm}<{\centering}| m{2cm}<{\centering} | m{0.8cm}<{\centering} |}' + '\n'
                                        + r'\arrayrulecolor{深蓝}\hline 药物 & 检测位点 & 检测结果 & 单体型 & 多态性类型 & 结论 & 等级 \\\hline' + '\n')
                                    line = re.sub(r'^&', r'\multirow{' + str(
                                        temp_chemsite[drug][site]['clinal'] + 1) + '}{*}{' + drug + '}&', line)
                                line = re.sub(r'label.haptype.' + drug + '.clin.' + gene,
                                              self.endchesite[drug][hap][site][2], line)
                                init_len += temp_chemsite[drug][site]['clinal']
                            if line.find('label.genetype.' + drug + '.clin.' + site) != -1:
                                # r'label.genetype.' + drug + '.clin.' + site
                                if init_len > total:
                                    # print('outline +++++++++++++++++++++++++----------')
                                    init_len = 0
                                    total = 34
                                    OUT.write('\end{xltabular}\n')
                                    OUT.write(r'\newpage' + '\n')
                                    OUT.write(
                                        r'\begin{xltabular}{\textwidth}{| m{1.8cm}<{\centering} | m{3.5cm}<{\centering} | m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.7cm}<{\centering}| m{2cm}<{\centering} | m{0.8cm}<{\centering} |}' + '\n'
                                        + r'\arrayrulecolor{深蓝}\hline 药物 & 检测位点 & 检测结果 & 单体型 & 多态性类型 & 结论 & 等级 \\\hline' + '\n')
                                    line = re.sub(r'^&', r'\multirow{' + str(
                                        temp_chemsite[drug][site]['conclusion']) + '}{*}{' + drug + '}&', line)
                                line = line.replace(r'label.genetype.' + drug + '.clin.' + site,
                                                    self.endchesite[drug][hap][site][2])
                                init_len += temp_chemsite[drug][site]['clinal']
                
                OUT.write(line)
    
    def percent(self, lines):
        wild = 0.33
        homozygote = 0.67
        pattern = re.compile('[0|1]/[1|2]:([\d]+),([\d]+)')
        ref, alt = pattern.search(lines[-1]).groups()
        per = 0
        if float(ref) + float(alt) != 0:
            per = float(alt) / (float(alt) + float(ref))
        if per < wild:
            result = lines[3] + lines[3]
        elif per > homozygote:
            result = lines[4] + lines[4]
        else:
            result = lines[3] + lines[4]
        print(result)
        return result


def main():
    parser = ArgumentParser()
    parser.add_argument('--chemhap', action='store', dest='chem', help='chemhaplotype.xls',
                        default='/mnt/archive/home/blood/Pipline/data/chemNew/chemhaplotype_0522_excel.xls')
    parser.add_argument('--drug', action='store', dest='drug', help='chemsite_drug.xls',
                        default='/mnt/archive/home/blood/Pipline/data/chemNew/chemsite_drug.xlsx')
    parser.add_argument('--anno', action='store', dest='anno', help='annofile multimutation', required=True)
    parser.add_argument('--C1', action='store', dest='tex', help='C1_infos.tex', required=True)
    parser.add_argument('--outC1', action='store', dest='out', help='C1_infos_temp.tex')
    parser.add_argument('--check', action='store', dest='checksite', help='outfile is SLCO1B1 597C>T  CC')
    arg = parser.parse_args()
    chemhaplotype = CHEMdrug(arg.chem, arg.drug, arg.anno, arg.tex)
    chemhaplotype.readhaptype()
    chemhaplotype.readdrug()
    chemhaplotype.searchdrug()
    chemhaplotype.outcheck(arg.checksite)
    chemhaplotype.modifyC1(arg.out)


if __name__ == '__main__':
    main()


