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
    def __init__(self, chemhap, drug, multimutaion, c1info,externalrs,out):
        self.haplotype = chemhap
        self.drug = drug
        self.c1info = c1info
        self.annofile = multimutaion
        self.haplotypedict = {}
        self.drugdict = {}
        self.endchesite = {}
        self.out=out
        self.external=externalrs
        global rslist

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
        global rslist
        work = xlrd.open_workbook(self.drug)
        rslist=[]
        sheet = work.sheet_by_index(0)
        nrows = sheet.nrows
        for row in range(1, nrows):
            rowvalue = sheet.row_values(row)
            if sheet.cell(row, 4).ctype == 0:
                ishaplotype = 'N'
            else:
                ishaplotype = 'H'
            Dict.addtodict5(self.drugdict,rowvalue[9],rowvalue[0], rowvalue[1], rowvalue[2]+':'+rowvalue[3], ishaplotype, rowvalue[5])
            rslist.append(rowvalue[3])

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
                if lines[mark] == 'rs2032582' and re.search(r'^1/2:',lines[-1]):
                    Dict.addtodict(mutationrs, lines[mark], 'CT')
                elif lines[mark] =='rs8175347':
                    [gt,fre] = self.percent(lines)
                    if gt=='TT':
                        nlty='(TA)6/(TA)6'
                    elif gt=='TA':
                        nlty='(TA)6/(TA)7'
                    elif gt=='AA':
                        nlty = '(TA)7/(TA)7'
                    Dict.addtodict(mutationrs, lines[mark], nlty)
                elif lines[mark] in rslist:
                    [gt,fre] = self.percent(lines)
                    Dict.addtodict(mutationrs, lines[mark], gt)
                else:
                    pattern1=re.compile('ABL1:NM_005157:exon\d+:c\.\w+:p\.(\w+)')
                    pattern2=re.compile('IDH2:NM_002168:exon\d+:c\.\w+:p\.(\w+)')
                    for rsp in rslist:
                        if pattern1.search(lines[9]):
                            paa=pattern1.search(lines[9]).groups(1)
                            if rsp == paa:
                                [gt,fre] = self.percent(lines)
                                Dict.addtodict(mutationrs,rsp, fre)
                        elif pattern2.search(lines[9]):
                            paa=pattern2.search(lines[9]).group(1)
                            if rsp == paa:
                                [gt,fre] = self.percent(lines)
                                Dict.addtodict(mutationrs, rsp, fre)
                        else:
                            pass
                    Dict.addtodict(mutationrs,'rs8175347','TA)6/(TA)6')
                            
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
        extern=self.external()
        mutationrs.update(extern)
        for partone in self.drugdict.keys():
            #########如果按部分写####
            for drug in self.drugdict[partone].keys():
                for gene in self.drugdict[partone][drug].keys():
                    rscom = []
                    rssite = {}
                    for rsid_c in self.drugdict[partone][drug][gene].keys():
                        [rsid,caa]=str(rsid_c).split(':')[0]
                        genesite=gene+'('+caa+')'
                        for hap in self.drugdict[partone][drug][gene][rsid_c].keys():
                            genetype = ''
                            if hap == 'N':
                                if rsid in mutationrs and str(rsid).startswith('rs'):
                                    genetype = rsid + ':' + mutationrs[rsid]
                                elif rsid in mutationrs and not str(rsid).startswith('rs'):
                                    genetype = '阳性'
                                elif str(rsid).startswith('rs'):
                                    genetype = rsid + ':' + self.drugdict[partone][drug][gene][rsid_c][hap]
                                elif not str(rsid).startswith('rs'):
                                    genetype = '阴性'

                                if genetype in self.haplotypedict[drug][gene].keys():
                                    Dict.addtodict4(self.endchesite, drug,gene,genesite,hap,
                                                    [genetype.split(':')[1]] + self.haplotypedict[drug][gene][genetype])
                                elif genetype=='阳性':
                                    Dict.addtodict4(self.endchesite, drug, gene, genesite,hap,[mutationrs[rsid]] + self.haplotypedict[drug][gene][rsid])
                                elif genetype=='阴性':
                                    Dict.addtodict4(self.endchesite,drug,gene, genesite,hap,[genetype]+self.haplotypedict[drug][gene][rsid])
                                else:
                                    print('{} u do not think about'.format(genetype))
                                    
                            else:
                                if rsid in mutationrs:
                                    genetype = rsid + ':' + mutationrs[rsid]
                                    rssite[genetype] = [genesite, mutationrs[rsid]]
                                else:
                                    genetype = rsid + ':' + self.drugdict[partone][drug][gene][rsid_c][hap]
                                    rssite[genetype] = [genesite, self.drugdict[partone][drug][gene][rsid_c][hap]]
                                rscom.append(genetype)
                if rscom:
                    combiners = self.comiter(rscom)
                    mark = 0
                    for coms in combiners:
                        if coms in self.haplotypedict[drug][gene].keys():
                            mark += 1
                            for site in coms.split(','):
                                if site in rssite:
                                    Dict.addtodict4(self.endchesite, drug, gene,rssite[site][0],'H',[rssite[site][1]]+self.haplotypedict[drug][gene][coms])
                                else:
                                    print('this rssite is not in list ' + site)
                        else:
                            pass
                    if mark == 0:
                        for site in rssite:
                            Dict.addtodict4(self.endchesite,drug,gene,rssite[site][0], 'H',
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
                for key3 in self.endchesite[key1][key2].keys():
                        for key4 in self.endchesite[key1][key2][key4].keys():
                            if key4=='N':
                                [basic, hploy, clinal, conclusion] = self.endchesite[key1][key2][key3][key4]
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
                for gene in self.endchesite[drug].keys():
                    for site in self.endchesite[drug][gene].keys():
                        for hap in self.endchesite[drug][gene][site].keys():
                            OUT.write(gene+'\t'+site+'\t'+self.endchesite[drug][gene][site][hap][0]+'\n')



    def modifyC1(self):
        temp_chemsite=self.computlines()
        OUT=open(self.out,'w')
        for c1tex in ['chemotherapy.tex','drug_target.tex','Support_therapy.tex']:
            parttex=os.path.join(self.c1info,c1tex)
            with open(parttex,'r') as T:
                for line in T:
                    for drug in self.endchesite.keys():
                        for gene in self.endchesite[drug].keys():
                            for site in self.endchesite[drug][gene].keys():
                                for hap in self.endchesite[drug][gene][site].keys():
                                    if re.search(r'{numline}{*}' + '{' + drug + '}', line):
                                        line = re.sub('{numline}', temp_chemsite[drug], line)
                                    elif re.search(r'label.haptype.'+drug+'.'+gene,line):
                                        line=re.sub(r'label.haptype.'+drug+'.'+gene,self.endchesite[drug][gene][site]['H'][1],line)
                                    elif re.search(r'label.toxicity.'+drug+'.'+gene,line):
                                        line=re.sub(r'label.toxicity.'+drug+'.'+gene,self.endchesite[drug][gene][site]['H'][3],line)
                                    elif re.search(site,line) and re.search(r'XX',line):
                                        line=re.sub(r'XX',self.endchesite[drug][gene][site][hap][0],line)
                                    elif re.search(r'label.haptype.'+drug+'.clin.'+gene,line):
                                        line=re.sub(r'label.haptype.'+drug+'.clin.'+gene,self.endchesite[drug][gene][site]['H'][2],line)
                                    elif re.search(r'label.genetype.'+drug+'.clin.'+site,line):
                                        line=re.sub(r'label.genetype.'+drug+'.clin.'+site,line,self.endchesite[drug][gene][hap][2],line)
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
        per=('%.2f\%%' % (per * 100))
        return [result,per]


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


