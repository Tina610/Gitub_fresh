#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/4/27 下午6:40
# @Author  : liting
# @contact : tinating610@163.com
# @File    : report_excel.py
# @Software: PyCharm
import os
import sys
import glob
import re
import xlrd
import itertools
from datetime import datetime
import argparse

pathscript = os.path.realpath(__file__)
binp = os.path.dirname(pathscript)
parser=argparse.ArgumentParser()
parser.add_argument('--genedata',help='gene database',action='store',dest='gene',default='gene2',required=True)
parser.add_argument('--mutationdata',help='mutation database',dest='muta',action='store',default='classficationxlsx',required=True)
parser.add_argument('--leuresult',help='leu result',action='store',dest='leu',required=True)
parser.add_argument('--combin',help='gene combination lib',action='store',dest='combin',required=True)
parser.add_argument('--inlatex',help='latex indir',action='store',dest='inlatex',required=True)
parser.add_argument('--outlatex',help='latex outdir result/out/file/',action='store',dest='out',required=True)
parser.add_argument('--genelist',help='vary project information',dest='genelist',action='store',required=True)
parser.add_argument('--infoexcel',help='sample list',action='store',dest='info',required=True)
args = parser.parse_args()

def main():
    inputfile = glob.glob(args.leu + '/*leu')[0]
    sample=os.path.basename(args.out)
    if not 'HB15' in sample:
        os._exit(0)
    else:
        report=data_basic(inputfile,args.gene,args.muta,args.combin,args.inlatex,args.out,args.genelist,args.info)
        report.run()
class data_basic():
    def __init__(self,inputfile,genedata,mutationdata,combin,inlatex,outlatex,genelist,info):
        self.infile=inputfile
        self.genedata=genedata
        self.mutemp=mutationdata
        self.combinfile=combin
        self.inlatex=inlatex
        self.outlatex=outlatex
        self.genelist=genelist
        self.infoexcel=info
        self.gene={}
        self.mutationdata={}
        self.combin={}

    def adddict(self,dictself,key_a, key_b, val):
        if key_a in dictself:
            dictself[key_a].update({key_b: val})
        else:
            dictself.update({key_a: {key_b: val}})

    def readfile(self):
        work=xlrd.open_workbook(self.mutemp)
        sheet=work.sheet_by_index(0)
        mark=int(0)
        nl=len(sheet.row_values(0))
        for idx,tl in enumerate(sheet.row_values(0)):
            if '突变类号' in tl:
                mark=int(idx)
            else:
                pass
        for col in range(mark+1,nl,2):
            for row in range(1,sheet.nrows):
                mutid = sheet.row_values(row)[mark]
                if sheet.cell(row,col).ctype==0 and sheet.cell(row,col+1).ctype==0:continue
                if sheet.cell(row,col).ctype==0:
                    clin='N'
                else:
                    clin=sheet.cell(row,col).value
                if sheet.cell(row,col+1).ctype==0:
                    valid='N'
                else:
                    valid=sheet.cell(row,col+1).value
                values=[clin,valid]
                self.adddict(self.mutationdata,mutid,sheet.cell(0,col).value,values)####添加证据级别值为list


    def readgene(self):
        mutaionfile=xlrd.open_workbook(self.genedata)
        information=mutaionfile.sheet_by_index(0)
        for row in range(1,information.nrows):
            keys=information.cell(row,0).value
            type2=information.cell(row,1).value###old is 3
            value1=re.sub('\n|\r','',str(information.cell(row,2).value))
            if value1=='':value1='基因简介无'
            value2=re.sub('\n|\r','',str(information.cell(row,3).value))
            if value2=='':value2='临床价值无'
            #val2=r'\\临床价值：'+str(value2)
            value=[value1,value2]   ####这一步可以对value1和value2为空滤过，看需求，目前写明"无"
            self.adddict(self.gene,keys,type2,value)

    def readcombin(self):
        combination=xlrd.open_workbook(self.combinfile)
        sheet=combination.sheet_by_index(0)
        for row in range(1,sheet.nrows):
            self.combin[sheet.cell(row,0).value]=sheet.cell(row,1).value

    def dealcebpa(self,cebpa,type):
        if cebpa:
            if len(cebpa)==1:
                value=self.mutationdata['xyb1CEB002'][type]
                self.genemerge['CEBPA']['xyb1CEB001'][0][5]=value[1]
                self.genemerge['CEBPA']['xyb1CEB001'][0][6]=value[2]
            else:
                pass
        else:
            pass
    def computeresult(self,gene):
        zhengju=0
        nline1=0
        if gene in self.genemerge.keys():
            num_nm=len(list(self.genemerge[gene].keys()))
            for muid in self.genemerge[gene].keys():
                for i in self.genemerge[gene][muid]:
                    string1,nline1=self.fxstr(i[-1],31)
                    string12,nline2=self.fxstr(i[-2],4)
                    zhengju=zhengju+nline2

        return  zhengju,nline1

    def readsample(self):
        name=os.path.basename(self.infile)
        genelist = []
        genecom = []
        cebpa = []
        self.genemerge = {}
        self.clinal ={}
        combine = ''
        type = ''
        patient = ''
        if 'N.leu' in name:
            with open(self.infile) as F:
                title=F.readline()
                type=F.readline().strip()
                patient=F.readline().strip()
                return combine,[type,patient]
        else:
            with open(self.infile,'r') as F:#,encoding="latin-1"
                title=F.readline()
                pattern1 = re.compile(r'(?P<nmid>NM\\_\w+)\((?P<gene>\w+)\):(?P<aa>c..+?)\((?P<chaa>.*?)\)')
                pattern2 = re.compile(r'(?P<nmid>NM\\_\w+)\((?P<gene>\w+)\):(?P<aa>c..+?)')
                #pattern2 = re.compile(r'(?P<nmid>NM\\_\w+)\((?P<gene>\w+)\):?(?P<aa>c..+)?')
                for line in F:
                    lines=line.strip('\n').split('\t')
                    exon=lines[1]
                    if exon=='splicing':
                        hgvs=pattern2.search(lines[0])
                    else:
                        hgvs=pattern1.search(lines[0])
                    idmu = lines[4]
                    type=lines[3]
                    patient=lines[5]
                    project=lines[6]
                    genelist.append(idmu)
                    nmid=hgvs.group('nmid')
                    gene=hgvs.group('gene')
                    caa=hgvs.group('aa')
                    caa,nlc=self.fxstr(caa,10)
                    paa=hgvs.group('chaa') if hgvs.group('chaa') else 'NULL'
                    paa,nlp=self.fxstr(paa,10)
                    # hgvs=[hgvs[i:i+46] for i in range(0,len(hgvs),46)]
                    if gene=='CEBPA' and idmu =='xyb1CEB001':
                        hgvs=lines[0]
                        cebpa.append([gene,lines[1],hgvs,lines[2]])
                    if type in self.mutationdata[idmu]:
                        valid=self.mutationdata[idmu][type][1]
                        prognosis = self.mutationdata[idmu][type][0]
                    if not gene in self.genemerge.keys():
                        self.genemerge[gene]={idmu:[]}
                        self.genemerge[gene][idmu].append([nmid,exon,caa,paa,lines[2],valid,prognosis])
                    else:
                        if not idmu in self.genemerge[gene].keys():
                            self.genemerge[gene][idmu].append([nmid,exon,caa,paa,lines[2],valid,prognosis])
                        else:
                            self.genemerge[gene][idmu].append([nmid,exon,caa,paa,lines[2],valid,prognosis])
                    if not gene in self.clinal and gene in self.gene:
                        if type in self.gene[gene]:
                            self.clinal[gene]=r'\\'.join(self.gene[gene][type])
                        else:
                            self.clinal[gene]=r'\\'.join(['基因检测无','临床价值无'])
                    else:
                        pass

            for i in range(2,6):
                for com in itertools.permutations(genelist,i):
                    andcom='&'.join(com)
                    if andcom in self.combin:
                        genecom.append(andcom)
                    else:
                        pass
            if genecom:
                print(genecom)
                end=max(genecom,key=len)
                genecom=self.combin[end]
            else:
                pass
            self.dealcebpa(cebpa, type)
            return genecom,[type,patient]

    def zh_chinese(self,strs):
        zh_pattern = re.compile(u'[\u4e00-\u9fff]+')
        if re.findall(zh_pattern, strs):
            return True
        else:
            return False

    def fxstr(self,str,le):
        strings = [str[i:i + le] for i in range(0, len(str), le)]
        lines=len(strings)
        strings=r'\\'.join(strings)
        string='\makecell{'+strings+'}'
        return string,lines

    def convertstr(self,chrstr,ln):
        strfirst=''
        longs=0
        j=0
        shang=len(chrstr.encode('UTF-8'))//ln
        mo=len(chrstr.encode('UTF-8'))%ln
        print('test====='+chrstr)
        for i in chrstr:
            unstr=i.encode('UTF-8')
            longs+=len(unstr)
            if longs>ln :
                j+=1
                strfirst=strfirst+r'\\'+i
                print(strfirst)
                longs=0
            else:
                strfirst=strfirst+i
        return strfirst

    def C1info(self):
        self.readfile()
        self.readgene()
        self.readcombin()
        combine,patient_info=self.readsample()
        Leu=self.outlatex+'/LeuReport'
        os.makedirs(Leu,exist_ok=True)
        if not self.genemerge:
            negtemp=self.inlatex+'/NegativeLeuReport/*'
            os.system('cp -rf '+negtemp+' '+Leu)
        else:
            postemp=self.inlatex+'/PositiveLeuReport/*'
            os.system('cp -rf '+postemp+' '+Leu)
        self.sample=os.path.basename(self.outlatex)
        getc12='/home/blood/bin/Python/bin/python3 {bin}/sampleC12.py -info {info} -genelist {genelist} -sam {sam} -C1 {outlatex}/LeuReport/chapters/C1_infos.tex -cls {inlatex}/PositiveLeuReport/med-report.cls -outc1 {outlatex}/LeuReport/chapters/C1_temp.tex -outmed {outlatex}/LeuReport/med-report.cls'\
            .format(bin=binp,info=self.infoexcel,genelist=self.genelist,sam=self.sample,inlatex=self.inlatex,outlatex=self.outlatex)
        print(getc12)
        os.system(getc12)
        OUT=open(self.outlatex+'/LeuReport/chapters/C1_infos.tex','w')
        with open(self.outlatex+'/LeuReport/chapters/C1_temp.tex','r') as F:
            lines=F.readlines()
            if '检测内容范围内，未检测到有临床意义突变' in lines[22]:
                os.system('cp {road}/LeuReport/chapters/C1_temp.tex {road}/LeuReport/chapters/C1_infos.tex'.format(road=self.outlatex))
            else:
                for line in lines:
                    if 'textcolor{海洋绿}{核苷酸改变}' in line:
                        OUT.write(line)
                        n=0
                        genes = list(sorted(self.genemerge.keys()))
                        firstgene = genes[0]
                        genes.pop(0)
                        # zhengju, yuhou = self.computeresult(firstgene)
                        # print(firstgene,zhengju,yuhou)
                        total = 8
                        for idmu in self.genemerge[firstgene]:
                            n+=1
                            samid=0
                            if len(self.genemerge[firstgene][idmu]) >1:
                                OUT.write(r'\multirow{2}{*}{\textbf{' + firstgene + r'}}')
                                for i in self.genemerge[firstgene][idmu]:
                                    samid+=1
                                    result=''
                                    zhengju,nline=self.fxstr(i[-2],4)
                                    print(zhengju,nline)
                                    for j in range(0,len(self.genemerge[firstgene][idmu][0])-2):
                                        result=result+r'&\textbf{'+i[j]+'}'
                                    if samid==2 and total==8:
                                        OUT.write(result + r'&\textbf{' + i[j + 1] + r'}\\\hline'+'\n')
                                        OUT.write(r'\end{xltabular}' + '\n')
                                        OUT.write(r'\newpage'+'\n')
                                        OUT.write(r'\arrayrulecolor{海洋绿}'+'\n')
                                        OUT.write(
                                            r'\begin{xltabular}{\textwidth}{| m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.5cm}<{\centering} | m{2cm}<{\centering} |X<{\centering}  | m{1cm}<{\centering} | m{1.5cm}<{\centering} |}\arrayrulecolor{海洋绿}\whline' + '\n')
                                        OUT.write(r'\rowcolor{浅灰色}' + '\n')
                                        OUT.write(
                                            r'\textcolor{海洋绿}{\normalsize{\textbf{\setlength{\baselineskip}{11pt}\makecell{突变\\基因}}}}&\normalsize{\textbf{\textcolor{海洋绿}{转录本ID}}}&\normalsize{\textbf{\textcolor{海洋绿}{突变位置}}}&\normalsize{\textbf{\textcolor{海洋绿}{核苷酸改变}}}&\normalsize{\textbf{\textcolor{海洋绿}{氨基酸改变}}}&\normalsize{\textbf{\setlength{\baselineskip}{11pt}\textcolor{海洋绿}{突变频率}}}&\normalsize{\setlength{\baselineskip}{11pt}\textbf{\textcolor{海洋绿}{\makecell{证据\\等级}}}} \\\whline' + '\n')
                                        result = ''
                                        total=37
                                    elif nline>2 and samid==1:
                                        print('this is n==1')
                                        OUT.write(result + r'&\textbf{' + i[j + 1] + r'}\\\hline')
                                        OUT.write(r'\end{xltabular}' + '\n')
                                        OUT.write(r'\newpage' + '\n')
                                        OUT.write(r'\arrayrulecolor{海洋绿}')
                                        OUT.write(
                                            r'\begin{xltabular}{\textwidth}{| m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.5cm}<{\centering} | m{2cm}<{\centering} |X<{\centering}  | m{1cm}<{\centering} | m{1.5cm}<{\centering} |}\arrayrulecolor{海洋绿}\whline' + '\n')
                                        OUT.write(r'\rowcolor{浅灰色}' + '\n')
                                        OUT.write(
                                            r'\textcolor{海洋绿}{\normalsize{\textbf{\setlength{\baselineskip}{11pt}\makecell{突变\\基因}}}}&\normalsize{\textbf{\textcolor{海洋绿}{转录本ID}}}&\normalsize{\textbf{\textcolor{海洋绿}{突变位置}}}&\normalsize{\textbf{\textcolor{海洋绿}{核苷酸改变}}}&\normalsize{\textbf{\textcolor{海洋绿}{氨基酸改变}}}&\normalsize{\textbf{\setlength{\baselineskip}{11pt}\textcolor{海洋绿}{突变频率}}}&\normalsize{\setlength{\baselineskip}{11pt}\textbf{\textcolor{海洋绿}{\makecell{证据\\等级}}}} \\\whline' + '\n')
                                        result=''
                                        total=37
                                    else:
                                        OUT.write(result + r'&\textbf{' + i[j + 1] + r'}\\\cline{2-7}')
                                OUT.write(r'&\multicolumn{6}{ p{12cm} |}{' + self.genemerge[firstgene][idmu][0][6] + '}')
                                if len(self.genemerge[firstgene].keys()) == n:
                                    OUT.write(r'\\\shline' + '\n')
                                else:
                                    OUT.write(r'\\\cline{2 - 7}' + '\n')

                            else:
                                zhengju, nline = self.fxstr(self.genemerge[firstgene][idmu][0][-2], 4)
                                yuhou,nline2=  self.fxstr(self.genemerge[firstgene][idmu][0][-1], 31)
                                print(nline,nline2,n)
                                if n==2 and total==8:
                                    OUT.write(r'\end{xltabular}' + '\n')
                                    OUT.write(r'\newpage' + '\n')
                                    OUT.write(r'\arrayrulecolor{海洋绿}')
                                    OUT.write(
                                        r'\begin{xltabular}{\textwidth}{| m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.5cm}<{\centering} | m{2cm}<{\centering} |X<{\centering}  | m{1cm}<{\centering} | m{1.5cm}<{\centering} |}\arrayrulecolor{海洋绿}\whline' + '\n')
                                    OUT.write(r'\rowcolor{浅灰色}' + '\n')
                                    OUT.write(
                                        r'\textcolor{海洋绿}{\normalsize{\textbf{\setlength{\baselineskip}{11pt}\makecell{突变\\基因}}}}&\normalsize{\textbf{\textcolor{海洋绿}{转录本ID}}}&\normalsize{\textbf{\textcolor{海洋绿}{突变位置}}}&\normalsize{\textbf{\textcolor{海洋绿}{核苷酸改变}}}&\normalsize{\textbf{\textcolor{海洋绿}{氨基酸改变}}}&\normalsize{\textbf{\setlength{\baselineskip}{11pt}\textcolor{海洋绿}{突变频率}}}&\normalsize{\setlength{\baselineskip}{11pt}\textbf{\textcolor{海洋绿}{\makecell{证据\\等级}}}} \\\whline' + '\n')
                                OUT.write(r'\multirow{2}{*}{\textbf{' + firstgene + r'}}')
                                result = ''
                                for i in range(0, len(self.genemerge[firstgene][idmu][0]) - 2):
                                    result = result + r'&\textbf{' + self.genemerge[firstgene][idmu][0][i] + '}'
                                OUT.write(result + r'&\textbf{' + self.genemerge[firstgene][idmu][0][i + 1])
                                if n==1 and nline>3 or nline2>3:
                                    OUT.write(r'}\\\hline')
                                    OUT.write(r'\end{xltabular}' + '\n')
                                    OUT.write(r'\newpage' + '\n')
                                    OUT.write(r'\arrayrulecolor{海洋绿}')
                                    OUT.write(
                                        r'\begin{xltabular}{\textwidth}{| m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.5cm}<{\centering} | m{2cm}<{\centering} |X<{\centering}  | m{1cm}<{\centering} | m{1.5cm}<{\centering} |}\arrayrulecolor{海洋绿}\whline' + '\n')
                                    OUT.write(r'\rowcolor{浅灰色}' + '\n')
                                    OUT.write(
                                        r'\textcolor{海洋绿}{\normalsize{\textbf{\setlength{\baselineskip}{11pt}\makecell{突变\\基因}}}}&\normalsize{\textbf{\textcolor{海洋绿}{转录本ID}}}&\normalsize{\textbf{\textcolor{海洋绿}{突变位置}}}&\normalsize{\textbf{\textcolor{海洋绿}{核苷酸改变}}}&\normalsize{\textbf{\textcolor{海洋绿}{氨基酸改变}}}&\normalsize{\textbf{\setlength{\baselineskip}{11pt}\textcolor{海洋绿}{突变频率}}}&\normalsize{\setlength{\baselineskip}{11pt}\textbf{\textcolor{海洋绿}{\makecell{证据\\等级}}}} \\\whline' + '\n')
                                    total=37
                                else:
                                    OUT.write(r'}\\\cline{2-7}')
                                OUT.write(r'&\multicolumn{6}{ p{12cm} |}{' + self.genemerge[firstgene][idmu][0][
                                    6] + '}')
                                if len(self.genemerge[firstgene].keys()) == n:
                                    OUT.write(r'\\\shline' + '\n')
                                else:
                                    OUT.write(r'\\\cline{2 - 7}' + '\n')

                        for gene in genes:
                            n=0
                            init=0
                            zhengju, yuhou = self.computeresult(gene)
                            print('this {} has {} {}'.format(gene,zhengju,yuhou))
                            if total==8:
                                OUT.write(r'\end{xltabular}' + '\n')
                                OUT.write(r'\newpage' + '\n')
                                OUT.write(r'\arrayrulecolor{海洋绿}')
                                OUT.write(
                                    r'\begin{xltabular}{\textwidth}{| m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.5cm}<{\centering} | m{2cm}<{\centering} |X<{\centering}  | m{1cm}<{\centering} | m{1.5cm}<{\centering} |}\arrayrulecolor{海洋绿}\whline' + '\n')
                                OUT.write(r'\rowcolor{浅灰色}' + '\n')
                                OUT.write(
                                    r'\textcolor{海洋绿}{\normalsize{\textbf{\setlength{\baselineskip}{11pt}\makecell{突变\\基因}}}}&\normalsize{\textbf{\textcolor{海洋绿}{转录本ID}}}&\normalsize{\textbf{\textcolor{海洋绿}{突变位置}}}&\normalsize{\textbf{\textcolor{海洋绿}{核苷酸改变}}}&\normalsize{\textbf{\textcolor{海洋绿}{氨基酸改变}}}&\normalsize{\textbf{\setlength{\baselineskip}{11pt}\textcolor{海洋绿}{突变频率}}}&\normalsize{\setlength{\baselineskip}{11pt}\textbf{\textcolor{海洋绿}{\makecell{证据\\等级}}}} \\\whline' + '\n')
                                total = 40
                            else:
                                total=37

                            if init>total:
                                OUT.write(r'\end{xltabular}' + '\n')
                                OUT.write(r'\newpage' + '\n')
                                OUT.write(r'\arrayrulecolor{海洋绿}')
                                OUT.write(
                                    r'\begin{xltabular}{\textwidth}{| m{1cm}<{\centering} | m{2cm}<{\centering} | m{1.5cm}<{\centering} | m{2cm}<{\centering} |X<{\centering}  | m{1cm}<{\centering} | m{1.5cm}<{\centering} |}\arrayrulecolor{海洋绿}\whline' + '\n')
                                OUT.write(r'\rowcolor{浅灰色}' + '\n')
                                OUT.write(
                                    r'\textcolor{海洋绿}{\normalsize{\textbf{\setlength{\baselineskip}{11pt}\makecell{突变\\基因}}}}&\normalsize{\textbf{\textcolor{海洋绿}{转录本ID}}}&\normalsize{\textbf{\textcolor{海洋绿}{突变位置}}}&\normalsize{\textbf{\textcolor{海洋绿}{核苷酸改变}}}&\normalsize{\textbf{\textcolor{海洋绿}{氨基酸改变}}}&\normalsize{\textbf{\setlength{\baselineskip}{11pt}\textcolor{海洋绿}{突变频率}}}&\normalsize{\setlength{\baselineskip}{11pt}\textbf{\textcolor{海洋绿}{\makecell{证据\\等级}}}} \\\whline' + '\n')
                                init=0
                            else:
                                init=init+zhengju+yuhou

                            for idmu in self.genemerge[gene]:
                                n+=1
                                if len(self.genemerge[gene][idmu]) >1:
                                    OUT.write(r'\multirow{2}{*}{\textbf{' + gene + r'}}')
                                    for i in self.genemerge[gene][idmu]:
                                        result=''
                                        for j in range(0,len(self.genemerge[gene][idmu][0])-2):
                                            result=result+r'&\textbf{'+i[j]+'}'
                                        OUT.write(result+r'&\textbf{'+i[j+1]+ r'}\\\cline{2-7}')
                                    OUT.write(r'&\multicolumn{6}{ p{12cm} |}{' +self.genemerge[gene][idmu][0][6]+'}')
                                    if len(self.genemerge[gene].keys())==n:
                                        OUT.write(r'\\\shline'+'\n')
                                    else:
                                        OUT.write(r'\\\cline{2 - 7}'+'\n')
                                else:
                                    OUT.write(r'\multirow{2}{*}{\textbf{' + gene + r'}}')
                                    result = ''
                                    for i in range(0,len(self.genemerge[gene][idmu][0])-2):
                                        result=result+r'&\textbf{'+self.genemerge[gene][idmu][0][i]+'}'
                                    OUT.write(result+r'&\textbf{'+self.genemerge[gene][idmu][0][i+1]+ r'}\\\cline{2-7}')
                                    OUT.write(r'&\multicolumn{6}{ p{12cm} |}{' + self.genemerge[gene][idmu][0][
                                            6] + '}')
                                    if len(self.genemerge[gene].keys()) == n:
                                        OUT.write(r'\\\hline' + '\n')
                                    else:
                                        OUT.write(r'\\\cline{2 - 7}' + '\n')

                    elif combine and re.search(r'\end{xltabular}',line):
                        OUT.write(line)
                        OUT.write(r'\item \normalsize{\textbf{'+'注:'+combine+'}}\n')
                    else:
                        OUT.write(line)
                os.system('rm {road}/LeuReport/chapters/C1_temp.tex'.format(road=self.outlatex))
            return patient_info


    def c2info(self):
        OUT2 = open(self.outlatex + '/LeuReport/chapters/C2_mutations_temp.tex', 'w')
        with open(self.outlatex + '/LeuReport/chapters/C2_mutations.tex', 'r') as F2:
            for line in F2:
                if re.findall('{genelabel}', line):
                    for key in sorted(self.clinal):
                        OUT2.write(r'\item {}\\'.format(key))
                        OUT2.write('\n')
                        OUT2.write(self.clinal[key] + r'\\' + '\n')
                else:
                    OUT2.write(line)
        os.system('mv {road}/LeuReport/chapters/C2_mutations_temp.tex {road}/LeuReport/chapters/C2_mutations.tex'.format(
            road=self.outlatex))


    def c3info(self):
        patient_info=self.C1info()
        dirtypes = os.listdir(self.outlatex + '/LeuReport/genelist')
        print(dirtypes)
        typeclass = 'C3_genelist_{}.tex'.format(patient_info[0])
        if typeclass in dirtypes:
            dirtype = os.path.join(self.outlatex + '/LeuReport/genelist', typeclass)
            print(dirtype)
            os.system(
                'mv {typegene} {road}/LeuReport/chapters/C3_genelist.tex'.format(typegene=dirtype, road=self.outlatex))
        else:
            raise Exception('do not have detect message')
        return patient_info[1]
    def run(self):
        samplename=self.c3info()
        self.c2info()
        today=datetime.today()
        curdate=today.strftime("%Y-%m-%d")
        reportname='{da}_{name}_{sample}-血液病'.format(da=curdate,name=samplename,sample=self.sample)
        os.system('cp {outlatex}/LeuReport/main.tex {outlatex}/LeuReport/{samid}.tex'.format(inlatex=self.inlatex,outlatex=self.outlatex,samid=reportname))
        os.system('cd {outlatex}/LeuReport/ && /home/blood/Software/textlive/2017/bin/x86_64-linux/xelatex -file-line-error -interaction=nonstopmode {outlatex}/LeuReport/{samid}.tex >2 /dev/null '.format(outlatex=self.outlatex,samid=reportname))
        os.system('cp {outlatex}/LeuReport/{samid}.pdf {outlatex}/'.format(outlatex=self.outlatex,samid=reportname))


if __name__ == '__main__':
    # temp=database('/Users/liting/Desktop/新报告模板/HB15JXEF00128-Positive_V1version.txt','/Users/liting/Desktop/新报告模板/金橡基因库模板.xls',
    #               '/Users/liting/Desktop/新报告模板/金橡突变库模板.xlsx','/Users/liting/Desktop/新报告模板/outresult.txt')
    # temp.run()
    main()
