#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/5/14 下午2:34
# @Author  : liting
# @contact : tinating610@163.com
# @File    : chemmodual.py
# @Software: PyCharm
import os, sys
import re
import Dict
import xlrd
import itertools
from argparse import ArgumentParser



class CHEMdrug():
    def __init__(self, chemhap, drug, multimutaion, c1info,externalrs):
        self.haplotype = chemhap
        self.drug = drug
        self.c1info = c1info
        self.annofile = multimutaion
        self.haplotypedict = {}
        self.drugdict = {}
        self.endchesite = {}
        self.external=externalrs

    def readhaptype(self):
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

    def readdrug(self):
        work = xlrd.open_workbook(self.drug)
        threeclass={}
        sheet = work.sheet_by_index(0)
        nrows = sheet.nrows
        for row in range(1, nrows):
            rowvalue = sheet.row_values(row)
            if sheet.cell(row, 4).ctype == 0:
                ishaplotype = 'N'
            else:
                ishaplotype = 'H'
            sitegene=rowvalue[1]+'('+rowvalue[2]+')'
            Dict.addtodict5(self.drugdict,rowvalue[9],rowvalue[0], sitegene, rowvalue[3], ishaplotype, rowvalue[5])

    def readannofile(self):
        mutationrs = {}
        listrs = self.rs2list()#######读取结果
        with open(self.annofile, 'r') as F:
            mark = None
            title = F.readline().strip().split('\t')
            for idx, value in enumerate(title):
                if 'snp13' in value:
                    mark = idx

            for line in F:
                lines = line.strip('\n').split('\t')
                if lines[mark] in listrs:
                    if lines[mark] == 'rs2032582' and re.search(r'^1/2:',lines[-1]):
                        Dict.addtodict(mutationrs, lines[mark], 'CT')
                    else:
                        gt = self.percent(lines)
                        Dict.addtodict(mutationrs, lines[mark], gt)

        return mutationrs
    
    def readexternal(self):
        extern={}
        work=xlrd.open_workbook(self.external)
        sheet=work.sheet_by_index(0)
        for row in range(1,sheet.nrows):
            rowvalue=sheet.row_values(row)
            if not rowvalue[1] in extern:
                extern[rowvalue[1]]=rowvalue[2]
            else:
                print('This rs'+rowvalue[1]+'has appears two')
                
        return extern
    

    def searchdrug(self):
        mutationrs = self.readannofile()
        for partone in ['part1','part2','part3']:
            #########如果按部分写####
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
                                        print('strange genetype: ' + hap+' '+rsid+' '+gene+drug)
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
        for key1 in self.endchesite.keys():
            drug_lines = 0
            for key2 in self.endchesite[key1].keys():
                if key2 == 'N':
                    for key3 in self.endchesite[key1][key2].keys():
                        [basic, hploy, clinal, conclusion] = self.endchesite[key1][key2][key3]
                        n_clinal = self.numlin(clinal, 35)
                        n_conclusion = self.numlin(conclusion, 3)
                        drug_lines = drug_lines + n_clinal * 0.7 + n_conclusion
                else:
                    key3list = list(self.endchesite[key1][key2].keys())
                    [basic, hploy, clinal, conclusion] = self.endchesite[key1][key2][key3list[0]]
                    n_clinal = self.numlin(clinal, 35)
                    drug_lines = drug_lines + len(key3list) + n_clinal * 0.7
            temp_chemsit[key1] = int(drug_lines)
        return temp_chemsit

    def outcheck(self,checkout):
        with open(checkout,'w') as OUT:
            for drug in self.endchesite:
                for hap in self.endchesite[drug].keys():
                    for site in self.endchesite[drug][hap].keys():
                        gene,cloc = re.search('(\w+)\((..+)\)', site).groups()
                        OUT.write(gene+'\t'+cloc+'\t'+self.endchesite[drug][hap][site][0]+'\n')



    def modifyC1(self,out):
        temp_chemsite = self.computlines()
        with open(out, 'w') as C1:
            F = open(self.c1info, 'r')
            for line in F:
                C1.write(line)
            F.close()

            if '氟达拉滨' in self.endchesite.keys():
                for key2 in self.endchesite['氟达拉滨'].keys():
                    for key3 in self.endchesite['氟达拉滨'][key2].keys():
                        values = self.endchesite['氟达拉滨'][key2][key3]
                        gene = re.search('(\w+)\(', key3).group(1)
                        [mutiltype, score] = self.multiclass['氟达拉滨'][gene]
                        rowfirst = temp_chemsite['氟达拉滨']
                        C1.write('\multirow{' + str(rowfirst) + '}{*}{氟达拉滨} & ' + key3 + '& ' + values[0] + ' & ' + values[1] + ' & ' + mutiltype + ' & ' + values[3] + ' & ' + score + r'\\\cline{2-7}& \multicolumn{6}{ p{13cm} |}{' + values[2] + r'}\\ \hline'+'\n')
                        C1.write(r'\end{xltabular}'+'\n')
                        C1.write(r'\newpage' + '\n')
            C1.write(
                r'\begin{xltabular}{\textwidth}{| m{1.8cm}<{\centering} | m{3.5cm}<{\centering} | m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.7cm}<{\centering}| m{2cm}<{\centering} | m{0.8cm}<{\centering} |}' + '\n'
                + r'\arrayrulecolor{深蓝}\hline 药物 & 检测位点 & 检测结果 & 单体型 & 多态性类型 & 结论 & 等级 \\\hline' + '\n')
            temp_chemsite.pop('氟达拉滨')
            self.endchesite.pop('氟达拉滨')
            init_len = 0
            while (len(temp_chemsite.keys()) > 0):
                newpage = False
                totalrow = 21
                for drug, nlines in sorted(temp_chemsite.items(), key=lambda item: item[1], reverse=True):
                    if init_len + int(nlines) > totalrow:
                        continue
                    else:
                        newpage = True
                    n = 0
                    if 'N' in self.endchesite[drug].keys() and 'H' in self.endchesite[drug].keys():
                        for hapy in sorted(self.endchesite[drug].keys()):
                            for site in self.endchesite[drug][hapy].keys():
                                n += 1
                                siterow = len(self.endchesite[drug][hapy].keys())
                                info = self.endchesite[drug][hapy][site]
                                gene = re.search('(\w+)\(', site).group(1)
                                if site in self.multiclass:
                                    [mutiltype, score] = self.multiclass[drug][site]
                                else:
                                    [mutiltype, score] = self.multiclass[drug][gene]

                                if hapy == 'H':
                                    if n == 1:
                                        C1.write('\multirow{' + str(nlines) + '}{*}{' + drug + '} & ' + site + ' & ' + info[0] + '& \multirow{' + str(siterow) + '}{*}{' + info[1] + '} & \multirow{' + str(siterow) + '}{2cm}{' + mutiltype + '} & \multirow{' + str(siterow) + '}{2cm}{' + info[3] + '} & \multirow{' + str(siterow) + '}{*}{' + score + r'} \\ \cline{2-3}' + '\n')
                                    elif n < siterow:
                                        C1.write('& ' + site + '&' + info[0] + '&&&&' + r'\\\cline{2-3}' + '\n')
                                    else:
                                        C1.write('& ' + site + '&' + info[0] + '&&&&' + r'\\\cline{2-7}' + '\n')
                                        C1.write('& \multicolumn{6}{ p{13cm} |}{' + info[2] + r'}\\\cline{2-7}' + '\n')

                                else:
                                    C1.write(r'& ' + site + '&' + info[0] + '& ' + info[1] + '& ' + mutiltype + '& ' + info[3] + '& ' + score + r'\\\cline{2-7}' + '\n')
                                    C1.write(r'& \multicolumn{6}{ p{13cm} |}{' + info[2] + r'}\\\hline' + '\n')
                        self.endchesite.pop(drug)
                        temp_chemsite.pop(drug)
                        init_len += nlines

                    elif ['N'] == list(self.endchesite[drug].keys()):
                        for site in self.endchesite[drug]['N']:
                            n += 1
                            siterow = len(self.endchesite[drug]['N'].keys())
                            info = self.endchesite[drug]['N'][site]
                            gene = re.search('(\w+)\(', site).group(1)
                            if site in self.multiclass:
                                [mutiltype, score] = self.multiclass[drug][site]
                            else:
                                [mutiltype, score] = self.multiclass[drug][gene]
                            if n == 1:
                                C1.write(
                                    r'\multirow{' + str(nlines) + '}{*}{' + drug + '} & ' + site + '& ' + info[0] + ' & ' +info[1] + '&' + mutiltype + '&' + info[3] + '&' + score + r'\\ \cline{2-7}' + '\n')
                                C1.write(r'& \multicolumn{6}{ p{13cm} |}{' + info[2])
                                if n==siterow:
                                    C1.write(r'}\\ \hline' + '\n')
                                else:
                                    C1.write(r'}\\ \cline{2-7}' + '\n')

                            else:
                                C1.write('& ' + site + '& ' + info[0] + '& ' + info[1] + '& ' + mutiltype + '& ' + info[3] + '& ' + score + r'\\\cline{2-7}' + '\n')
                                C1.write(r'& \multicolumn{6}{ p{13cm} |}{' + info[2])
                                if n == siterow:
                                    C1.write(r'}\\ \hline' + '\n')
                                else:
                                    C1.write(r'}\\\cline{2-7}' + '\n')
                        self.endchesite.pop(drug)
                        temp_chemsite.pop(drug)
                        init_len += nlines

                    elif ['H'] == list(self.endchesite[drug].keys()):
                        for site in self.endchesite[drug]['H']:
                            n += 1
                            siterow = len(self.endchesite[drug]['H'].keys())
                            info = self.endchesite[drug]['H'][site]
                            gene = re.search('(\w+)\(', site).group(1)
                            if site in self.multiclass:
                                [mutiltype, score] = self.multiclass[drug][site]
                            else:
                                [mutiltype, score] = self.multiclass[drug][gene]
                            if n == 1:
                                C1.write('\multirow{' + str(nlines) + '}{*}{' + drug + '} & ' + site + ' & ' + info[0] + '& \multirow{' + str(siterow) + '}{*}{' + info[1] + '} & \multirow{' + str(siterow) + '}{2cm}{' + mutiltype +'} & \multirow{' + str(siterow) + '}{2cm}{' + info[3] + '} & \multirow{' + str(siterow) + '}{*}{' + score + r'} \\ \cline{2-3}' + '\n')
                            elif n < siterow:
                                C1.write('& ' + site + '&' + info[0] + '&&&&' + r'\\\cline{2-3}' + '\n')
                            else:
                                C1.write('& ' + site + '&' + info[0] + '&&&&' + r'\\\cline{2-7}' + '\n')
                                C1.write('& \multicolumn{6}{ p{13cm} |}{' + info[2] + r'}\\\hline' + '\n')
                        self.endchesite.pop(drug)
                        temp_chemsite.pop(drug)
                        init_len += nlines
                if not newpage:
                    C1.write(r'\end{xltabular}' + '\n')
                    C1.write(r'\newpage' + '\n')
                    C1.write(
                        r'\begin{xltabular}{\textwidth}{| m{1.8cm}<{\centering} | m{3.5cm}<{\centering} | m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.7cm}<{\centering}| m{2cm}<{\centering} | m{0.8cm}<{\centering} |}' + '\n'+r'\arrayrulecolor{深蓝}\hline 药物 & 检测位点 & 检测结果 & 单体型 & 多态性类型 & 结论 & 等级 \\\hline' + '\n')
                    init_len = 0
                else:
                    pass
            C1.write(r'\end{xltabular}' + '\n')
            C1.write(r'\newpage' + '\n')


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
    parser.add_argument('--chemhap', action='store', dest='chem', help='chemhaplotype.xls', default='/home/liting/Leukemia/PiplineReport/data/chemNew/chemhaplotype.xls')
    parser.add_argument('--drug', action='store', dest='drug', help='chemsite_drug.xls', default='/home/liting/Leukemia/PiplineReport/data/chemNew/chemsite_drug.xlsx')
    parser.add_argument('--anno', action='store', dest='anno', help='annofile multimutation', required=True)
    parser.add_argument('--C1', action='store', dest='tex', help='C1_infos.tex', required=True)
    parser.add_argument('--outC1',action='store',dest='out',help='C1_infos_temp.tex')
    parser.add_argument('--external',action='store',dest='extern',help='haved the result',required=True)
    parser.add_argument('--check',action='store',dest='checksite',help='outfile is SLCO1B1 597C>T  CC')
    arg = parser.parse_args()
    chemhaplotype = CHEMdrug(arg.chem,arg.drug,arg.anno,arg.tex)
    chemhaplotype.readhaptype()
    chemhaplotype.readdrug()
    chemhaplotype.searchdrug()
    chemhaplotype.outcheck(arg.checksite)
    chemhaplotype.modifyC1(arg.out)

if __name__ == '__main__':
    main()


