#!/usr/bin/perl -w

## MISC_getHelp.pl -ifile "C:/Perl/site/lib/funclib.pl" -ofile "e:/funclib.help.txt"
#  mean                     
#  SS                       
#  cor                      
#  getPermute               
#  fileMode                 
#  ncolFile                 
#  isFileEmptyOrNotExist      说明: 文件为空或不存在时返回1
#  getFileName                说明: 返回文件名
#  getExt                     说明: 返回扩展名,只认最后一个.
#  getDir                     说明: 返回路径
#  writeLog                 
#  currTime                   说明: 返回当前时间字符串,格式2005-03-18 08:56:38
#  NONPSSM                    说明: PSSM值不存在时的值
#  getTmpPath                 说明: 仅限于本机的临时路径. bar=1末端有/,0无/
#  value2key                  说明: 由value得到key名 hash{key}=value
#  hashlength                 说明: 得到哈希表长度
#  insertLine2File            说明: 在文件中第几行插入,插入后该行行号为row,首1
#  trim                     
#  replaceStr                 说明: 将str内的某些字符改为新字符串,若没找到,则不替换
#  shellSort                  说明: 若为asc则升序,其它字符都为降序
#  sortInt                    说明: 若为asc则升序,其它字符都为降序
#  getFileNames               说明: 返回文件名包含路径,输入路径可含/或不含,匹配扩展名要用如"\.txt"或 txt$;
#  saveMtx2File               说明: 保存矩阵到文件,默认覆盖,如果成功,返回1.
#  loadFile2Mtx               说明: 加载文件到矩阵中,跳过skipn行.
#  loadFile2String            说明: 加载文件到字符串中,跳过skipn行.
#  sample                     说明: 将ref拆分成nbin份,对每份在source中采样.
#  flipMatrix                 说明: 转置矩阵,都是$类型的矩阵.
#  mtx2str                  
#  findOverlaps               说明: 对两个矩阵,查找重叠情况, 参考R中的IRanges.findOverlaps()
#  countOverlaps              说明: 对两个矩阵,计算qry在sbj中的个数, 参考R中的IRanges.countOverlaps()
#  getIntervals             
#  fillGaps                 
#  seqFormat                  说明: 判断序列文件是fa或fq格式,只通过第1行的标题是>还是@来判断格式
#  isBad                      说明: qc,判断序列是否含 10%的N 或 QT%的ATCG..
#  remove12                 
#  splitRaw                 
#  getUniqByCols            
#  splitFileByCols          
#  findTail                   说明: 从seq的from(不包含from位置,首1)开始连续的找A或T,允许首个A/T离from的4个位置以内.
#  trimseq                    说明: 传longseq的引用,$isRC=1需要反转互补,0不需要
#  reverseAndComplement       说明: $r=1 reverse $c=1 complement
#  trimseqFromTo              说明: 传longseq的引用,$isRC=1需要反转互补,0不需要
#  grpSame                  
#  grpByPos_I               
#  grpByPos                 
#  isIP                       说明:不判断边界,默认位点左右各有9,10nt (-10~-1[PA],1~10)
#  formatPatOutput            说明:格式化patronus的输出,并替换原文件
#  getNt                      说明: 根据背景概率(整数的),产生联子
#  convertChrToIdx            说明: atcgATCG-->下标
#  convertIdxToChr            说明:下标-->ATCG
#  getKgramId                 说明:根据kgram,得到下标
#  genOneKGram                说明:产生1个空的k联子，idx表示下标
#  genKgrams                  说明:产生空的k联子，withvalue表示是否右边产生＝0
#  cntKgramsByK               说明:统计全部k-gram在seqfile中出现次数,$from,$to<1则为统计整条序列
#  cntKgrams                  说明:统计给定grams在seqfile中出现次数,$from,$to<1则为统计整条序列
#  cntEachPosByK              说明:
#  cntEachPosByGrams          说明:统计seqfile中从from~to的位置的各kgram出现次数,必须设置from,to
#  kcnt2pssm                  说明:kcnt阵转换为pssm阵
#  sortKcnt                   说明:降序排序kcnt矩阵,以所有列的和排序
#  sortPssm                   说明:降序排序pssm矩阵,以所有列的最大值排序
#  fas2tbl                    说明: 读入fas,输出2列matrix(title,seq)
#  mutateMotif                说明: 最后结果在@mutes中,每个值1个motif
#  execSql                    说明:执行sql语句
#  cloneTblIndex              说明: 在totbl中创建fromtbl的index
#  tblExists                  说明:判断table是否存在
#  connectDB                  说明:连接DB
#  getFldsIdx                 说明:返回表中指定字段的下标
#  getTblFlds                 说明:返回表中所有字段名
#  getFldValues               说明:取得某列的值,返回数组. ncol from 0.
#  loadFile2Tbl               说明:文件导入数据表
#  getTblFlds_lite            说明:返回表中所有字段名
#  getFldsIdx_lite            说明:返回表中指定字段的下标
#  file2LiteTbl               说明:文件导入Lite数据表
#  sql2file                 
#  getSqlFlds_Lite          
#  createPAtbl              
#  geneFromGff              



use strict;
use File::Basename;
use File::Copy;
use List::Util qw(shuffle);
use List::MoreUtils qw/ uniq /; 
use XML::Parser;
use XML::Simple;
use DBI;
use Time::Local;
use Text::NSP::Measures::2D::Fisher::twotailed;


#############################################################################
#公共函数库 require ("funclib.pl");
#############################################################################

#****************************************************************************
# 公共变量 
#****************************************************************************
#用于findOverlap中的type
our ($ANY,$WITHIN,$CONTAIN,$EQUAL,$OVP)=(1,2,3,4,5);

our $XMLDIR='E:/sys/code/XML_conf/';
our $TMPDIR='E:/projectData/';

##用our xx在脚本中声明。
our $PL_MAPGFF='E:/sys/code/PAT/PAT_mapGff.pl';
our $PL_ALTERPAs='E:/sys/code/PAT/PAT_alterPAs.pl';
our $PL_PA2PAC='E:/sys/code/PAT/PAT_PA2PAC.pl';
our $PL_DOAMB='E:/sys/code/UTIL/UTIL_DoAMB.pl';

our $RICEJPCOLS="dry_seed1:dry_seed2:dry_seed3;embryo1:embryo2;endosperm1:endosperm2:endosperm3;imbibed_seed1:imbibed_seed2:imbibed_seed3;shoot1:shoot2:shoot3;leaf_20days1:leaf_20days2:leaf_20days3;leaf_60days1:leaf_60days2:leaf_60days3;stem_60days1:stem_60days2:stem_60days3;root_5days1:root_5days2:root_5days3;root_60days1:root_60days2:root_60days3;husk1:husk2:husk3;anther1:anther2:anther3;mature_pollen1:mature_pollen2:mature_pollen3;pistil1:pistil2:pistil3";
our $INDICACOLS="dry_seed1:dry_seed2:dry_seed3;embryo1:embryo2:embryo3;endosperm1:endosperm2:endosperm3;imbibed_seed1:imbibed_seed2:imbibed_seed3;shoot1:shoot2:shoot3;leaf_20days1:leaf_20days2:leaf_20days3;leaf_60days1:leaf_60days2:leaf_60days3;stem_60days1:stem_60days2:stem_60days3;root_5days1:root_5days2:root_5days3;root_60days1:root_60days2:root_60days3;husk1:husk2:husk3;anther1:anther2:anther3;mature_pollen1:mature_pollen2:mature_pollen3;pistil1:pistil2:pistil3";
our $RICETWOCOLS="i_dry_seed1:i_dry_seed2:i_dry_seed3;i_embryo1:i_embryo2:i_embryo3;i_endosperm1:i_endosperm2:i_endosperm3;i_imbibed_seed1:i_imbibed_seed2:i_imbibed_seed3;i_shoot1:i_shoot2:i_shoot3;i_leaf_20days1:i_leaf_20days2:i_leaf_20days3;i_leaf_60days1:i_leaf_60days2:i_leaf_60days3;i_stem_60days1:i_stem_60days2:i_stem_60days3;i_root_5days1:i_root_5days2:i_root_5days3;i_root_60days1:i_root_60days2:i_root_60days3;i_husk1:i_husk2:i_husk3;i_anther1:i_anther2:i_anther3;i_mature_pollen1:i_mature_pollen2:i_mature_pollen3;i_pistil1:i_pistil2:i_pistil3;j_dry_seed1:j_dry_seed2:j_dry_seed3;j_embryo1:j_embryo2;j_endosperm1:j_endosperm2:j_endosperm3;j_imbibed_seed1:j_imbibed_seed2:j_imbibed_seed3;j_shoot1:j_shoot2:j_shoot3;j_leaf_20days1:j_leaf_20days2:j_leaf_20days3;j_leaf_60days1:j_leaf_60days2:j_leaf_60days3;j_stem_60days1:j_stem_60days2:j_stem_60days3;j_root_5days1:j_root_5days2:j_root_5days3;j_root_60days1:j_root_60days2:j_root_60days3;j_husk1:j_husk2:j_husk3;j_anther1:j_anther2:j_anther3;j_mature_pollen1:j_mature_pollen2:j_mature_pollen3;j_pistil1:j_pistil2:j_pistil3";
our $RICEJPCOLSj="j_dry_seed1:j_dry_seed2:j_dry_seed3;j_embryo1:j_embryo2;j_endosperm1:j_endosperm2:j_endosperm3;j_imbibed_seed1:j_imbibed_seed2:j_imbibed_seed3;j_shoot1:j_shoot2:j_shoot3;j_leaf_20days1:j_leaf_20days2:j_leaf_20days3;j_leaf_60days1:j_leaf_60days2:j_leaf_60days3;j_stem_60days1:j_stem_60days2:j_stem_60days3;j_root_5days1:j_root_5days2:j_root_5days3;j_root_60days1:j_root_60days2:j_root_60days3;j_husk1:j_husk2:j_husk3;j_anther1:j_anther2:j_anther3;j_mature_pollen1:j_mature_pollen2:j_mature_pollen3;j_pistil1:j_pistil2:j_pistil3";
our $INDICACOLSi="i_dry_seed1:i_dry_seed2:i_dry_seed3;i_embryo1:i_embryo2:i_embryo3;i_endosperm1:i_endosperm2:i_endosperm3;i_imbibed_seed1:i_imbibed_seed2:i_imbibed_seed3;i_shoot1:i_shoot2:i_shoot3;i_leaf_20days1:i_leaf_20days2:i_leaf_20days3;i_leaf_60days1:i_leaf_60days2:i_leaf_60days3;i_stem_60days1:i_stem_60days2:i_stem_60days3;i_root_5days1:i_root_5days2:i_root_5days3;i_root_60days1:i_root_60days2:i_root_60days3;i_husk1:i_husk2:i_husk3;i_anther1:i_anther2:i_anther3;i_mature_pollen1:i_mature_pollen2:i_mature_pollen3;i_pistil1:i_pistil2:i_pistil3";

our %convertor = (
    'TCA' => 'S', # Serine
    'TCC' => 'S', # Serine
    'TCG' => 'S', # Serine
    'TCT' => 'S', # Serine
    'TTC' => 'F', # Phenylalanine
    'TTT' => 'F', # Phenylalanine
    'TTA' => 'L', # Leucine
    'TTG' => 'L', # Leucine
    'TAC' => 'Y', # Tyrosine
    'TAT' => 'Y', # Tyrosine
    'TAA' => '_', # Stop
    'TAG' => '_', # Stop
    'TGC' => 'C', # Cysteine
    'TGT' => 'C', # Cysteine
    'TGA' => '_', # Stop
    'TGG' => 'W', # Tryptophan
    'CTA' => 'L', # Leucine
    'CTC' => 'L', # Leucine
    'CTG' => 'L', # Leucine
    'CTT' => 'L', # Leucine
    'CCA' => 'P', # Proline
    'CCC' => 'P', # Proline
    'CCG' => 'P', # Proline
    'CCT' => 'P', # Proline
    'CAC' => 'H', # Histidine
    'CAT' => 'H', # Histidine
    'CAA' => 'Q', # Glutamine
    'CAG' => 'Q', # Glutamine
    'CGA' => 'R', # Arginine
    'CGC' => 'R', # Arginine
    'CGG' => 'R', # Arginine
    'CGT' => 'R', # Arginine
    'ATA' => 'I', # Isoleucine
    'ATC' => 'I', # Isoleucine
    'ATT' => 'I', # Isoleucine
    'ATG' => 'M', # Methionine
    'ACA' => 'T', # Threonine
    'ACC' => 'T', # Threonine
    'ACG' => 'T', # Threonine
    'ACT' => 'T', # Threonine
    'AAC' => 'N', # Asparagine
    'AAT' => 'N', # Asparagine
    'AAA' => 'K', # Lysine
    'AAG' => 'K', # Lysine
    'AGC' => 'S', # Serine
    'AGT' => 'S', # Serine
    'AGA' => 'R', # Arginine
    'AGG' => 'R', # Arginine
    'GTA' => 'V', # Valine
    'GTC' => 'V', # Valine
    'GTG' => 'V', # Valine
    'GTT' => 'V', # Valine
    'GCA' => 'A', # Alanine
    'GCC' => 'A', # Alanine
    'GCG' => 'A', # Alanine
    'GCT' => 'A', # Alanine
    'GAC' => 'D', # Aspartic Acid
    'GAT' => 'D', # Aspartic Acid
    'GAA' => 'E', # Glutamic Acid
    'GAG' => 'E', # Glutamic Acid
    'GGA' => 'G', # Glycine
    'GGC' => 'G', # Glycine
    'GGG' => 'G', # Glycine
    'GGT' => 'G', # Glycine
    );


#****************************************************************************
# APA
#****************************************************************************
##由样本字段，判断是否是 RICE系，不是的话，直接返回原串
sub getOptionSmp {
  my $smp=shift;
  #RICEJPCOLS,INDICACOLS,RICETWOCOLS
  if ($smp=~/^RICEJPCOLS$/i) {
	$smp=$RICEJPCOLS;
  } elsif ($smp=~/^INDICACOLS$/i) {
	$smp=$INDICACOLS;
  } elsif ($smp=~/^RICETWOCOLS$/i) {
	$smp=$RICETWOCOLS;
  }elsif ($smp=~/^RICEJPCOLSj$/i) {
	$smp=$RICEJPCOLSj;
  }elsif ($smp=~/^INDICACOLSi$/i) {
	$smp=$INDICACOLSi;
  }
  return($smp);
}

#****************************************************************************
# Stat
#****************************************************************************
##计算log2
sub log2 {
  my $n = shift;
  return log($n)/log(2);
}

## 计算均值
## $m=mean(\@data)
sub mean {
my ($a)=@_;
return(0) if $#$a==-1;
my ($i,$sum)=(0,0);
for ($i=0;$i<=$#$a;$i++){
  $sum=$sum+$$a[$i];
}
return($sum/($#$a+1));
}

## 计算标准差
## $s=SS(\@data)
sub SS {
my ($a)=@_;
my ($i,$sum)=(0,0);
my $m=mean($a);
for ($i=0;$i<=$#$a;$i++){
  $sum=$sum+($$a[$i]-$m)*($$a[$i]-$m);
}
return(sqrt(($sum/$#$a)));
}

## 计算pearson相关性
## $c=cor(\@x,\@y);
## 如果x或y的标准差为0,则返回-2
sub cor {
my($x,$y)=@_;
die "not equal number for cor" if ($#$x!=$#$y or $#$x==-1);
my $mx=mean($x);
my $my=mean($y);
my $sx=SS($x);
my $sy=SS($y);
return(-2) if ($sx==0 or $sy==0);
my $sum=0;
for my $i(0..$#$x) {
  $sum+=($$x[$i]-$mx)*($$y[$i]-$my)/$sx/$sy;  
}
return ($sum/$#$x);
}

#############################################################################
## 给定矩阵，得到所有排列情况，输出矩阵
## 参考：http://topic.csdn.net/u/20110324/13/5c108144-8b05-47ac-8853-cfd97c7ed8a6.html?seed=1145643556&r=72355220
## $a=getPermute($mtx);
## 如$mtx=([1,2],[3,4,5])，则得到 ([1,3],[1,4],...[2,5])
#############################################################################
sub getPermute { 
  my $mtx=shift;
  my ($i,$m,$k);
  my $ret=[];
  for $i (0..$#$mtx){
    if($#$ret == -1 ){
	  for my $ii(0..$#{$mtx->[$i]}) {
		push(@$ret,[$mtx->[$i][$ii]]);
	  }        
       next;
    }
    my $temp =[];
    for $m(0..$#$ret){
	  for $k(0..$#{$mtx->[$i]}) {
		push(@$temp,[@{$ret->[$m]},$mtx->[$i][$k]]);
	  }
    }
    $ret = $temp;	
  }#i
return $ret;
}

#############################################################################
## Fisher检验
## $pval=fisherTest(n11,n12,n21,n22);
#--------------------------------------
#          word2   ~word2
#  word1    *n11      n12 | *n1p
# ~word1     n21      n22 |  n2p
#           --------------
#           *np1      np2   *npp
#############################################################################
sub fisherTest { 
  my ($n11,$n12,$n21,$n22)=@_;
	my $n1p=$n11+$n12;
	my $np1=$n11+$n21;
	my $npp=$n11+$n12+$n21+$n22;
    my $pvalue=calculateStatistic( n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
	return($pvalue);
}

#****************************************************************************
# Common
#****************************************************************************
#############################################################################
## $mode=fileMode($file)
## 取得文件的类型: win,unix,mac,none
## 参考 http://www.phpfans.net/article/htmls/201003/MjgyMTk1.html 
#############################################################################
sub fileMode {
    my $file = shift;    
    if ( -s $file and -T $file )  {
        open (TEST, "$file") or return undef;
        binmode(TEST);
        my $stream = "\0" x 1024;
        read ( TEST, $stream, 1024 );
        close TEST;
		return "win" if ( $stream =~ /\015\012/o );
		return "unix" if ( $stream =~ /[^\015]\012/o );
		return "mac" if ( $stream =~ /\015[^\012]/o );
		return "NONE";
    }  elsif ( ! -e $file )  {
        return "NONE";
    }
 }

#############################################################################
#mysep=getSep($aFileName)
#取得文件的分隔符，是\s或\t
#############################################################################
sub getSep {
  my($filename)=shift;
  open(FILE, $filename) or die "Could not open file '$filename' $!";
  my $line=<FILE>;
  close(FILE);
  chomp($line);
  my @seps=('\t','\s',';',',');
  for my $sep(@seps) {
	  my @elements = split ($sep, $line);
	  if (scalar(@elements)>1) {
		return($sep);
	  }
  }
  die "getSep: cannot find sep for $filename\n";
}

#############################################################################
## $nc=ncolFile($file)
## 取得文件的列数,只判断第1行
#############################################################################
sub ncolFile {
  my $file=shift;
  return(0) if !defined($file);
  ##确定列数
  my $ncol=0;
  open(ININ,"<$file") or die "Cannot read $file";
  my $line='';
  $line=<ININ>;
  $line=trim($line);
  if ($line ne '') {
    my @nn=split(/\t/,$line);
    $ncol=scalar(@nn);
  }
  close(ININ);
  return($ncol);
}

#############################################################################
#  isFileEmptyOrNotExist(afilename) 
#  useage: isFileEmptyOrNotExist('xx/xx.txt'):0/1       
#  说明: 文件为空或不存在时返回1
#  NOTE: 即使文件名中带空格 $f="c:/dir dir/h h.txt"也可以判断正确.
#############################################################################
sub isFileEmptyOrNotExist {
  my $f=shift;
  return(1) if (!(-e $f) or (-z $f));
  return(0);
}

#############################################################################
#  getFileName(afilename,noExt=0/1) 
#  useage: getFileName('xx/xx.txt'); getFileName('xx/xx.txt',1)        
#  说明: 返回文件名（不含路径）
#############################################################################
sub getFileName {
  my $f=shift;
  my $noExt=shift;
  if (!$noExt) {
	return basename($f);
  } else {
	$f=basename($f);
	my $a=rindex($f,'.');
	if ($a==-1) {
	  return $f;
	} else {
	  return substr($f,0,$a);
	}
  }  
}

#############################################################################
#  getExt(afilename) 
#  useage: getExt('xx/xx.txt')       
#  说明: 返回扩展名,只认最后一个.
#############################################################################
sub getExt {
  my $f=shift;
  my $a=rindex($f,'.');
	if ($a==-1) {
	  return '';
	} else {
	  return substr($f,$a+1,length($f)-$a);
	}
}

#############################################################################
#  getDir(afilename,nobar=0/1) 
#  useage: getDir('xx/xx.txt'); getFileName('xx/xx.txt',1)        
#  说明: 返回路径
#############################################################################
sub getDir {
  my $f=shift;
  my $nobar=shift;
  my ($name,$dir,$suffix) = fileparse($f);
  if (!$nobar) {
	if ($dir eq '.') {
	  $dir.='/';
	}
  } else {
    if ($dir ne '.' and substr($dir,length($dir)-1,1) eq '/') {
	  $dir=substr($dir,0,length($dir)-1);
    }
  }
  return $dir;
}

sub writeLog {
  my($lf,$str,$timer)=@_;   
  if ($timer) {
	$str=currTime().' '.$str;
  }
  print $lf $str;
}

#############################################################################
#  currTime
#  useage: currTime()       
#  说明: 返回当前时间字符串,格式2005-03-18 08:56:38
#  2011/5/15 若GMT＝1,则返回GMT的时间戳
#############################################################################
sub currTime {
  my $GMT=shift;
  if (defined($GMT) and $GMT==1) {
	return timegm(localtime);
  }
  my   @array=(localtime)[5,4,3,2,1,0];   
  $array[0]+=1900;   
  $array[1]+=1;   
  my   $nowtime=sprintf("%04u-%02u-%02u %02u:%02u:%02u",@array);   
  return($nowtime);
}

#############################################################################
#  NONPSSM:-99
#  useage: NONPSSM()       
#  说明: PSSM值不存在时的值
#############################################################################
sub NONPSSM {
  return -99;
}

#############################################################################
#   getTmpPath($bar):str
#  useage: $str=getTmpPath(1)       
#  说明: 仅限于本机的临时路径. bar=1末端有/,0无/
#############################################################################
sub getTmpPath {
  my($bar)=shift;
  my $mine=$TMPDIR;
  if (-e $mine) {
    return($mine.'/') if $bar;
	return($mine);
  } else {
	my $others='c:/PALAB';
	mkdir($others) if !(-e $others);  
	return($others.'/') if $bar;
	return($others);
  }
}

#############################################################################
#  str2hash($str,$sep)
#  useage: $hash=str2hash('x1=xx;x2=xx;',';')
#  说明: 由'x1=xx;x2=xx;'串，生成哈希变量 $hash{X1}=xx...
#  2015/7/10 key全用大写
#############################################################################
sub str2hash { 
  my $str=shift;
  my $sep=shift;
  $sep=';' if !defined($sep);
  my %hash=();
  my @items=split($sep,$str);
  for my $item(@items) {
	my @x=();
	if ($item=~/=/) { #xx=yy的形式
	  @x=split('=',$item);
	} else {
	  @x=split('\s\"',$item); #genebiotype "protein_coding" (mm10.gtf)
	}
    if (scalar(@x)==2) {
	  $x[1]=~s/\"//g;
	  #print scalar(@x)."..$x[0]...$x[1]\n";
	  $hash{trim(uc($x[0]))}=trim($x[1]);
    }
  }
  return %hash;
}


#############################################################################
#  value2key($value,%fmap)
#  useage: key=value2key('xx',%fmap)
#  说明: 由value得到key名 hash{key}=value
#############################################################################
sub value2key {
  my($value,%fmap)=@_;
  my($k);
  foreach $k(keys(%fmap)) {
	return($k) if $fmap{$k} eq $value;
  }
}

#############################################################################
#  l=hashlength(%fmap)
#  useage: l=hashlength(%fmap)       
#  说明: 得到哈希表长度
#############################################################################
sub hashlength {
  my %fmap=@_;
  return(scalar(keys(%fmap)));
}

#############################################################################
#  insertLine2File($file,$line,$rowNum)
#  useage: insertLine2File('xx.txt','aaaa',1);       
#  说明: 在文件中第几行插入,插入后该行行号为row,首1
#############################################################################
sub insertLine2File {
  my($file,$line,$rowNum)=@_;
  my($l);
  my($f2)=$file."x";
  return if ($rowNum<=0);
  open (IN,"<$file");
  open(OUT,">$f2");
  while ($l=<IN>) {
	if ($.==$rowNum) {
      print OUT "$line\n"
	}
	print OUT $l;
  }
  close(IN);
  close(OUT);
  unlink($file);
  rename($f2,$file);
}


#############################################################################
#  function trim delete the whitespace at the begining and end of the string#
#  useage: $var=trim($var);                                                 #
#############################################################################
sub trim($)
{
	my $string = shift;
	return '' if (!defined($string) or $string eq '');
	$string =~ s/^\s+//;
	$string =~ s/[\s\r\n]+$//;
	return $string;
}

#############################################################################
#line=readLastLine(file)
#读文件最后一行
#把指针移动到最后面,一个一个字节往前读,直到读到\N为止
#文件的末尾只能是含有1个空行，或者没有空行，否则返回是空行！
#返回的line不含回车
#############################################################################
sub readLastLine {
	my $file = shift or die;
	my $content = "";

	open (F, $file) or die $!;
	seek (F, 0, 2);   # set handler at the end of $file
	until ($content =~ m/\n(.*)\n?$/)	{
			my $string;
			if (seek (F, -1024, 1))  # backward 1024 bytes
			{
					my $n = read (F, $string, 1024) or die $!;
					$content = $string . $content;
					last if ($n < 1024);
					seek(F, -1024, 1);
			}
			else
			{
					my $len = tell F;
					seek (F, 0, 0) || die "see error at $file";
					read (F, $string, $len) or die $!;
					$content = $string . $content;
					last;
			}
	}
	close(F);
	if ($content =~ m/\n\n$/) 	{
		return("\n");
	}
	elsif ($content =~ m/\n?(.*)\n?$/)	{
		return("$1");
	}
	else	{
		return($content);
	}

}

#############################################################################
#  replaceStr($str,$substr,$to):$newstr
#  useage: $str=replaceStr($str,'.','.n')
#  说明: 将str内的某些字符改为新字符串,若没找到,则不替换
#############################################################################
sub replaceStr {
  my($str,$substr,$to)=@_;
  return $str if index($str,$substr)==-1;
  substr($str,index($str,$substr),length($substr))=$to;
  return $str;
}

#############################################################################
#  round($number):$integer
#  useage: round(2.4)=2
#  说明: 四舍五入
#############################################################################
sub round {
my($number) = shift;
return int($number + .5);
}

#############################################################################
#  shellSort($first,$last,$order,@list):@order
#  useage: @order=shellSort(0,$#list,'desc',@list);
#  说明: 若为asc则升序,其它字符都为降序
#############################################################################
sub shellSort {
  my($first,$last,$asc1,@list)=@_;	
  my($i,$j,$h,$to,$ninth,$tv,$cmp);
  $h=1;
  my(@order)=(0..$#list);
  my $asc=0;
  $asc=1 if (lc($asc1) eq 'asc');
  $ninth=int(($last-$first)/9);
  while ($h <= $ninth) {
    $h =($h * 3) + 1;
  }
    while ($h > 0){
       for $i ($first + $h..$last) {
        $tv=$list[$i];
		$to=$order[$i];
        $j = $i;        
        while ($j >= $first+$h and ( ($tv-$list[$j-$h]<0 and $asc==1) or ($tv-$list[$j-$h]>0 and $asc!=1) ) ){
          $list[$j] = $list[$j-$h];
		  $order[$j] = $order[$j-$h];
          $j-=$h;
        }
        $list[$j] = $tv;
		$order[$j] = $to;
       } #for
    $h=int($h/3);
  } #while
  return @order;
}

#############################################################################
#  sortInt($mtx,$intCol,order):$omtx
#  useage: $omtx =sortInt($mtx,0); $omtx =sortInt($mtx,1,'desc');
#  说明: 若为asc则升序,其它字符都为降序
#  intCol为需要排序的列
#############################################################################
sub sortInt {
	my ($in,$col,$order)=@_;
	$order='asc' if !$order;
	my $out;
	if ($order eq 'asc') {
	  @$out = sort{$a->[$col] <=> $b->[$col]} @{$in};
	} else {
      @$out = sort{$b->[$col] <=> $a->[$col]} @{$in};
	}	
	return $out;
}

#############################################################################
#  sortStr($mtx,$strCol,order):$omtx
#  useage: $omtx =sortStr($mtx,0); $omtx =sortStr($mtx,1,'desc');
#  说明: 若为asc则升序,其它字符都为降序
#  intCol为需要排序的列
#############################################################################
sub sortStr {
	my ($in,$col,$order)=@_;
	$order='asc' if !$order;
	my $out;
	if ($order eq 'asc') {
	  @$out = sort{$a->[$col] cmp $b->[$col]} @{$in};
	} else {
      @$out = sort{$b->[$col] cmp $a->[$col]} @{$in};
	}	
	return $out;
}

#############################################################################
#  getFileNames($dir,$pattern,$NOTMATCH):@files
#  useage: @files=getFileNames('.','\.pl'); @files=getFileNames('.','\.pl',1);      
#  说明: 返回文件名包含路径,输入路径可含/或不含,匹配扩展名要用如"\.txt"或 txt$; 
#       若$NOTMATCH非０,则表示不匹配模式才输出文件
#  若要匹配 xxtrain1.arff; train222.arff; 则用 pat="train.*\.arff"
#  若要匹配以train开头以.arff结尾的;则用 pat="^train.*\.arff$"
#  以arff结尾的 arff$
#############################################################################
sub getFileNames{
	my ($dir)=shift;
	my ($pattern)=shift;
	my ($not)=shift;
	#print "$dir,$pattern\n";
	$pattern='' if !$pattern;
	$dir.='/' if (substr($dir, length($dir)-1, 1) ne '/');
	my(@files,$file);
	local *DH;
	if (!opendir(DH, $dir)) {
	  warn "Cannot opendir $dir: $! $^E";
	  next;
	}
	foreach (readdir(DH)) {
	  if ($_ eq '.' || $_ eq '..') {
	  next;
	  }
	  $file = $_; 
	  #print "$file\n";
	  if(-f $dir.$file) {
	    push(@files, $dir.$file) if ((!$not and $file=~/$pattern/i) or ($not and $file!~/$pattern/i));
	  }
	}
	closedir(DH);
	return @files;
}

#############################################################################
#  moveFiles($dir,$pattern,$destDir):0/1
#  useage: $success=moveFiles('.','\.pl','E:/")  
#  说明: 移动文件  
#  若要匹配 xxtrain1.arff; train222.arff; 则用 pat="train.*\.arff"
#  若要匹配以train开头以.arff结尾的;则用 pat="^train.*\.arff$"
#  返回 1 成功； 0 失败
#############################################################################
sub moveFiles{
	my ($dir)=shift;
	my ($pattern)=shift;
	my ($destDir)=shift;
	return 0 if !defined($destDir);
	return 0 if !(-e $dir);
	return 0 if !(-e $destDir);
	my(@files)=getFileNames($dir,$pattern);
	foreach my $file (@files) {
	  move("$file","$destDir")||warn "could not copy files :$!" ;
	}
	return 1;
}

#############################################################################
#  saveMtx2File($mtx,$file,$append,$blankstr,$idxfrom,$idxto):sucess
#  useage: if saveMtx2File($mtx,$file) {..} or saveMtx2File($mtx,$file,1)
#  saveMtx2File($mtx,$file,1,'\N',100,1000)
#  说明: 保存矩阵到文件,默认覆盖,如果成功,返回1.
#  saveMtx2File($mtx,$file,1,'\N') --对空值填充\N
#  真tmd危险,这里的FH,居然会引用到其它pl的FH去.最好取个BT点的句柄名
#  2012-09-10 增加$idxfrom,$idxto选项，可以选择输出mtx的某些行
#############################################################################
sub saveMtx2File {
  my($mtx,$file,$append,$blankstr,$idxfrom,$idxto)=@_;
  return(0) if $#$mtx==-1;
  if (!$append) {
	open(FH_saveMtx2File, ">$file") || return 0;  
  } else{
    open(FH_saveMtx2File, ">>$file") || return 0;  
  }  
  $idxfrom=0 if !defined($idxfrom) or $idxfrom<0;
  $idxto=$#$mtx if !defined($idxto) or $idxto>$#$mtx;
  die "saveMtx2File: idxfrom=$idxfrom; idxto=$idxto" if $idxfrom>$idxto;
  my($i,$j,$cols);
  for $i ( $idxfrom .. $idxto) {
	$cols=$#{$mtx->[$i]};
    for $j ( 0 .. ($cols-1) ) { 
	 if (!defined($mtx->[$i][$j])) {
	  print FH_saveMtx2File $blankstr."\t";
	 } else {
	  print FH_saveMtx2File "$mtx->[$i][$j]\t"; 	
	 }
  } 	
  if (!defined($mtx->[$i][$cols])) {
    print FH_saveMtx2File $blankstr."\n"; 
  } else {
    print FH_saveMtx2File "$mtx->[$i][$cols]\n";	
  }
} 
  close(FH_saveMtx2File);  
  return 1;
}


#############################################################################
#  loadFile2Mtx($file,$sep,$skipn):$mtx
#  useage: $mtx=loadFile2Mtx($file,"\t",1)
#  说明: 加载文件到矩阵中,跳过skipn行.
#  加快速度
#  2016/11/21 不得不再重新用trim -- 因为chomp就算放两个也没法消除window下的换行符，奇怪了b
#############################################################################
sub loadFile2Mtx {
  my($file,$sep,$n)=@_;
  my($mtx);
  $n=0 if !defined($n);
  open my $infh, '<', $file or return; 
  if ($n>0) {
	 while ($n>0) {
		<$infh>;
		$n--;
	 }
  }
  if (!defined($sep) or !$sep or $sep eq "\t" or $sep eq '\t') {
	  while ( my $line=<$infh> ) {
			$line=trim($line);
			next unless length($line); 
			push(@$mtx,[ split(/\t/,$line)]);
	  }
	  return($mtx);
  }

  while ( my $line=<$infh> ) {
			$line=trim($line);
		next unless length($line); 
		push(@$mtx,[ split($sep,$line)]);#用 $sep比直接'\t慢一倍！
  }
  return($mtx);
}


#############################################################################
#  loadFile2String($file,$skipn):$str
#  useage: $str=loadFile2String($file,1)
#  说明: 加载文件到字符串中,跳过skipn行.
#############################################################################
sub loadFile2String {
  my($file,$n)=@_;
  my($str);
  open(LOADTO, "<$file") || return;  
  my($line);
  while ($line=<LOADTO>) {
    next if $.<=$n;
	$str.=$line;
   }
  close(LOADTO);  
  return $str;
}

#############################################################################
#  sample($ref,$source,$nbin):@index
#  useage: @idx=sample(\@ref,\@source,5);
#  说明: 将ref拆分成nbin份,对每份在source中采样. 
#        返回的idx个数与ref相同,若某组无法得到相同个数,则返回空集.
#        !!输入的ref和source必须是已经排好序的. 输出的index是source的下标
#############################################################################
sub sample {
  my ($ref,$source,$nbin)=@_; #已经排好序了
  my ($i,$j,$min,$max,@grp,$id);
  my(@idx)=();
  my($pi)=0;
  my ($nr)=scalar(@$ref);
  my ($ns)=scalar(@$source);
  return(@idx) if ($ns<$nr or $nbin==0);
  
  #每组个数
  my ($ng)=int($nr/$nbin);
  my ($ngs)=$ng;
  my ($res)=0;
  if ($nr%$nbin) {	
	$res=$nr%$nbin;
	$nbin++;
  }

  #对每组的ref,取最小值和最大值,对source遍历得到在最小值最大值之间的值,如果个数>要求,则再随机抽取
  for $i(0..($nbin-1)) {
	@grp=();
    #最后一组
    if ($res and $i==$nbin-1) {
	  $ng=$res;
    }
	$min=$max=$$ref[$i*$ngs];
	for $j(1..($ng-1)) {
	  #最大最小值
	  $id=$i*$ngs+$j;
	  if ($min>$$ref[$id]) {
		$min=$$ref[$id];
	  }elsif ($max<$$ref[$id]) {
		$max=$$ref[$id];
	  }
	}
	#print "$ng,$min,$max\n";
	#遍历
    while ($pi<$ns) {
	  if ($$source[$pi]<$min) {
		$pi++;
	  } elsif ($$source[$pi]<=$max and $$source[$pi]>=$min) {
		push(@grp,$pi); #保存index
		$pi++;
      }elsif ($$source[$pi]>$max) {
        last;
      }	  
    }
	#判断是否满足个数
	#print "max=$max,min=$min,grp=@grp\n";
	if (scalar(@grp)<$ng) { #若少于,则退出
	  @idx=();
	  return(@idx);
	}elsif (scalar(@grp)>$ng) { #多于,则随机取
	  @grp=shuffle(@grp);
	  push(@idx,@grp[0..($ng-1)]);
	}else {
	  push(@idx,@grp);
	}	
  } #for i
  return(@idx);
}


#############################################################################
#  flipMatrix($matrix):$flipmatrix
#  useage: $fliparray=flipMatrix($array);
#  说明: 转置矩阵,都是$类型的矩阵.
#############################################################################
sub flipMatrix {
  my($array)=shift;
  my($x);
  my @fliparray =
	map { $x = $_;
	[ map { $$array[$_][$x] } 0 .. $#$array ];
	} 0 .. $#{$$array[0]};
  my $fliparray2=\@fliparray;
  return($fliparray2);
}

#############################################################################
#  mtx2str($mtx,$sep):$str
#  useage: $str=mtx2str($mtx,"\t");
#  矩阵转换为字符串,主要用于调试用
#############################################################################
sub mtx2str {
  my($mtx)=shift;
  my $sep=shift;
  $sep="\t" if (!defined($sep));
  my($i);
  my $str='';
  for $i ( 0 .. $#{$mtx} ) {
	$str.=join($sep,@{$mtx->[$i]})."\n";
  }
 return $str;
}

#############################################################################
#  findOverlaps($qry,$sbj,$qsCol,$qeCol,$ssCol,$seCol,$leftMargin,$rightMargin,$select,$drop,$type,$minOverlap,$outputType):$mtxIdx[qi,si]
#  注意: 是按qsCol和ssCol都排序后的!!!!
#  说明: 对两个矩阵,查找重叠情况, 参考R中的IRanges.findOverlaps()
#  $qry,$sbj: 待比较的2个矩阵,将qry比对到sbj中,返回的IDX矩阵行数与qry相同; 比如qry是PA,sbj是GFF
#  $qsCol,$qeCol,$ssCol,$seCol: start/end的列号(首0),若start=end,则表示是一个点,如PA
#  $leftMargin/rightMargin: 允许将sbj首尾扩展一定长度,若<0,则表示缩小一定长度
#  $select=('all','first') 当一个qry可以比对到多个sbj时(如一个PA比对到多个GFF)
#    first/last(未实现)时,选择的是第1个sbj或最后1个sbj,返回IDX矩阵列数为2; all则输出所有sbj的ID,这时返回IDX矩阵的列数不同.
#  $drop=1/0(default): 是否丢弃比对不到的qry; 1则丢弃. 比如将PA比对到GFF,若不需要比对不到的PA,则drop=1
#  $type=('any','within','contain','equal') any则当qry与sbj有重叠时都算,within是qry在sbj中,equal是完全相同
#  $minOverlap=1: 最小重叠长度,对type都有效.
#  $outputType=0(default)/1: 若1,则输出overlapType,主要用于当type=any的时候,输出不同的type
#  返回$mtxIdx[qi,si]: 矩阵行数取决于drop=1/0,矩阵列数取决于select=all/first 以及 outputType(=1,则第2列是type列(首1))

#  Usages:
#  1) 将PA比对到GFF中,只考虑比对到的PA,以及只取首个GFF: select=first; drop=1; type=any
#  PAT_mapGff.pl 中的mode=1
#  idx=findOverlaps($PA,$GFF,3,3,3,4,0,0,'first',1,'any',1)

#  2) 将PA (start/end)比对到GFF中,且记录覆盖类型
#  PAT_mapGff.pl 中的mode=2
#  idx=findOverlaps($PA,$GFF,3,4,3,4,0,0,'first',1,'any',1,1)
#############################################################################
sub findOverlaps {
  my($qry,$sbj,$qsCol,$qeCol,$ssCol,$seCol,$leftMargin,$rightMargin,$select,$drop,$type,$minOverlap,$outputType)=@_;
  $select=lc($select);
  $type=lc($type);
  die "Error select=$select" if ($select ne 'all' and $select ne 'first'); # and $select ne 'last'
  die "Error type=$type" if ($type ne 'any' and $type ne 'within' and $type ne 'equal' and $type ne 'contain');
  die "Error drop=$drop" if ($drop!=0 and $drop!=1);
  $minOverlap=1 if $minOverlap<=0;
  $type=$ANY if $type eq 'any';
  $type=$WITHIN if $type eq 'within';
  $type=$CONTAIN if $type eq 'contain';
  $type=$EQUAL if $type eq 'equal';
  my ($ALL,$FIRST,$LAST)=(1,2,3);
  $select=$ALL if $select eq 'all';
  $select=$FIRST if $select eq 'first';
  #$select=$LAST if $select eq 'last';
  
  my $mtxIdx;
  my $qLast=$#$qry;
  my $sLast=$#$sbj;
  my ($qi,$si)=(0,0);
  my ($qs,$qe,$ss,$se);
  my ($have,$ovpType,$oType)=(0,-1); #(未找到,不重合)
  my $ovpLen=0;
  my $lastSi=0;
  while ($qi<=$qLast) {
	$ovpType=-1;
	$have=0;
    $qs=$qry->[$qi][$qsCol];
	$qe=$qry->[$qi][$qeCol];
    $si=$lastSi; #当select=all时,会跳过si,这里需要根据上一个qi的首个si作为lastSi

	while ($si<=$sLast) {
	  $ss=$sbj->[$si][$ssCol];
	  $se=$sbj->[$si][$seCol];
	  $ss-=$leftMargin;
	  $se+=$rightMargin;
	  if ($se<$qs) {
		$si++;
        while ($si>$sLast and !$drop and $qi<=$qLast) { ## 2011/10/13 非drop时,需要补上si完了但qi没完的情况!!
          push(@$mtxIdx,[$qi]);
	      $qi++;
	      next;
        }
		next;
	  }
	  if ($ss>$qs) {
		if (!$have and !$drop) { #drop=1,即使没找着,也输出
		  push(@$mtxIdx,[$qi]);
		}
        last;
	  }
	  $lastSi=$si if $ovpType==-1;
	  #取得qry和sbj的匹配情况
      if ($qs==$ss and $qe==$se) {
		$ovpType=$EQUAL;
		$ovpLen=$qe-$qs+1;
      } elsif ($qs>=$ss and $qe<=$se) {
		$ovpType=$WITHIN;
		$ovpLen=$qe-$qs+1;
	  } elsif ($qs<=$ss and $qe>=$se) {
		$ovpType=$CONTAIN;
		$ovpLen=$se-$ss+1;
	  } elsif ($qs>=$ss and $qs<=$se) {
        $ovpType=$OVP;
		$ovpLen=$se-$qs+1;		
	  } elsif ($qe>=$ss and $qe<=$se) {
	    $ovpType=$OVP;
		$ovpLen=$qe-$ss+1;
	  }
	  $oType=$ovpType;
	  $ovpType=$ANY if $ovpType!=-1 and $type==$ANY;
      
	  #输出
	  if ($type==$ovpType and $ovpLen>=$minOverlap) {		
		if ($select==$FIRST) {
		  $have=1;
		  if ($outputType) {
			push(@$mtxIdx,[$qi,$oType,$si]);
		  } else {
			push(@$mtxIdx,[$qi,$si]);
		  }		  
		  last;
		}elsif ($select==$ALL) {
          if (!$have) {
			if ($outputType) {
			  push(@$mtxIdx,[$qi,$oType,$si]);
			} else {
              push(@$mtxIdx,[$qi,$si]);
		    }
			$have=1;
          } else {
			push(@{$mtxIdx->[$#$mtxIdx]},$si); 
		  }
		  $si++;
		  next;
	    }		
	  } else {
		push(@$mtxIdx,[$qi]) if !$drop; #虽然有重合,但type不对
		$si++;
		next;
	  }
	} #while si
  
   $qi++;
  } #while qi

  if (!$drop and $#$qry!=$#$mtxIdx) {
  	 saveMtx2File($qry,"findOverlaps.qry.txt");
	 saveMtx2File($sbj,"findOverlaps.sbj.txt");
	 saveMtx2File($mtxIdx,"findOverlaps.mtxIdx.txt");
     die "Error in findOverlaps(): drop=0 but nrow(qry) not eq nrow(mtxIdx); save input to findOverlaps.qry/sbj/mtxIdx.txt";
  }

  return $mtxIdx;
}

#############################################################################
#  countOverlaps($qry,$sbj,$qsCol,$qeCol,$ssCol,$seCol,$sumCols,$leftMargin,$rightMargin,$select,$drop,$type,$minOverlap):$mtxIdx
#  注意: 是按qsCol和ssCol都排序后的!!!!
#  说明: 对两个矩阵,计算qry在sbj中的个数, 参考R中的IRanges.countOverlaps()
#  $qry,$sbj: 待比较的2个矩阵,将qry比对到sbj中,返回的IDX矩阵行数与sbj相同; 比如qry是PA,sbj是GFF
#  $qsCol,$qeCol,$ssCol,$seCol: start/end的列号(首0),若start=end,则表示是一个点,如PA
#  $sumCols: 用于sum的列,若有多列则如 3:4:5, 比如 sum(leaf),sum(seed)
#  $leftMargin/rightMargin: 允许将sbj首尾扩展一定长度,若<0,则表示缩小一定长度
#  $select=('all','first') 当一个qry可以比对到多个sbj时(如一个PA比对到多个GFF)
#    first时,一个qry只计数一次;all时,一个qry可重复计数.比如PA比对到gene,若PA1同时在gene1和gene2中,则PA1可以计数2次(如果select=all)
#  $drop=1/0(default): 是否丢弃无qry的sbj; 1则丢弃. 比如将PA比对到GFF,若不需要比对不到的GFF,则drop=1
#  $type=('any','within','contain','equal') any则当qry与sbj有重叠时都算,within是qry在sbj中,equal是完全相同
#  $minOverlap=1: 最小重叠长度,对type都有效.
#  返回$mtxIdx[si,qryCnt,sumCols]: 矩阵行数与sbj同行，但取决于drop=1/0,矩阵列数取决于sumCols,注意第2列为qryCnt,如比对到gene的PA个数

#  Usages:
#  1) 将PA比对到GFF中,一个PA只能属于一个GFF,计数leaf/seed(2,3列),去除比对不到的GFF: select=first; drop=1; type=any; sumCols=2:3
#  idx[si,qryCnt,sumCols]=countOverlaps($PA,$GFF,3,3,3,4,'2:3',0,0,'first',1,'any',1)
#############################################################################
sub countOverlaps {
  my($qry,$sbj,$qsCol,$qeCol,$ssCol,$seCol,$sumCols,$leftMargin,$rightMargin,$select,$drop,$type,$minOverlap)=@_;
  $select=lc($select);
  $type=lc($type);
  die "Error select=$select" if ($select ne 'all' and $select ne 'first'); # and $select ne 'last'
  die "Error type=$type" if ($type ne 'any' and $type ne 'within' and $type ne 'equal' and $type ne 'contain');
  die "Error drop=$drop" if ($drop!=0 and $drop!=1);
  $minOverlap=1 if $minOverlap<=0;
  $type=$ANY if $type eq 'any';
  $type=$WITHIN if $type eq 'within';
  $type=$CONTAIN if $type eq 'contain';
  $type=$EQUAL if $type eq 'equal';
  my ($ALL,$FIRST,$LAST)=(1,2,3);
  $select=$ALL if $select eq 'all';
  $select=$FIRST if $select eq 'first';
  #$select=$LAST if $select eq 'last';
  my @sc=split(/:/,$sumCols); #要求和的列
  my $mtxIdx;
  my $qLast=$#$qry;
  my $sLast=$#$sbj;
  my ($qi,$si)=(0,0);
  my ($qs,$qe,$ss,$se);
  my ($have,$ovpType,$cnt)=(0,-1); #(未找到,不重合)
  my $ovpLen=0;
  my $lastQi=0;
  my @noQry=@sc;
  my @sums=@sc;
  for my $i(0..$#noQry) {
    $noQry[$i]=0;
	$sums[$i]=0;
  }
  while ($si<=$sLast) {
	$ovpType=-1;
	$have=$cnt=0;
	$ss=$sbj->[$si][$ssCol];
	$se=$sbj->[$si][$seCol];
	$ss-=$leftMargin;
	$se+=$rightMargin;
    if ($select == $ALL) { #如果qry可重复,则每次从前一个sbj开始的qry开始滑动
	  $qi=$lastQi;
    }	
	@sums=@noQry;	
	while ($qi<=$qLast) {
      $qs=$qry->[$qi][$qsCol];
	  $qe=$qry->[$qi][$qeCol];
	  if ($qe<$ss) {
		$qi++;
		$lastQi=$qi;
		next;
	  }

	  if ($qs>$se) {
        last;
	  }
	  #取得qry和sbj的匹配情况
      if ($qs==$ss and $qe==$se) {
		$ovpType=$EQUAL;
		$ovpLen=$qe-$qs+1;
      } elsif ($qs>=$ss and $qe<=$se) {
		$ovpType=$WITHIN;
		$ovpLen=$qe-$qs+1;
	  } elsif ($qs<=$ss and $qe>=$se) {
		$ovpType=$CONTAIN;
		$ovpLen=$se-$ss+1;
	  } elsif ($qs>=$ss and $qs<=$se) {
        $ovpType=$OVP;
		$ovpLen=$se-$qs+1;		
	  } elsif ($qe>=$ss and $qe<=$se) {
	    $ovpType=$OVP;
		$ovpLen=$qe-$ss+1;
	  }
	  $ovpType=$ANY if $ovpType!=-1 and $type==$ANY;      
	  #输出
	  if ($type==$ovpType and $ovpLen>=$minOverlap) {
		$have=1;
		$cnt++;
		for my $i(0..$#sc) {
		  $sums[$i]+=$qry->[$qi][$sc[$i]];
		}
	  }
      $qi++;
	} #while qi
	if (!$have and !$drop) { #drop=1,即使没找着,也输出
	  push(@$mtxIdx,[$si,0,@noQry]);
	} elsif ($have) {
      push(@$mtxIdx,[$si,$cnt,@sums]);
	}
   $si++;
  } #while si
  return $mtxIdx;
}


#############################################################################
#  getIntervals($mtx/$ary,$col,$more):$mtx
#  注意: 是排序后的!!!!
#  可以得到一串的start-end
#  ary: 已排序的数组,里面含有相同的元素;
#  mtx,col: 已按col列排序的矩阵
#  如果未定义col,则认为$ary=\@x形式
#  如果有定义col,则认为$mtx形式
#  输出$mtx: 每个唯一值所在的区间 (startIdx,endIdx)
#  2012-09-09 改进：$col允许多列，以1:2的形式提供; 增加选项$more=0(default)/1，如果是1,则输出的第3列为如chr1,-
#  2016/3/23 改进：增加 $more=N1:N2的形式('0:10')，如果提供，则只对 N1~N2行进行区间划分，获得idx（比如允许对某个gene范围内的tr进行区间划分，而不是对所有基因）
#Usage: 
#1) 输入数组
#my @ary=(1,1,1,2,2,3);
#my $ret2=getIntervals(\@ary); #不能指定col,否则出错

#2) 输入$型矩阵
#my $ret=getIntervals($mtx,0); #必须指定col

#3) 输入多个列
#my $ret=getIntervals($mtx,'0:1:3');

#4) 只对给定行范围进行区间划分
#my $ret=getIntervals($mtx,5,0:10); #只对0-10行的数据进行划分 
#############################################################################
sub getIntervals {
	my($mtx)=shift;
	my $col=shift;
	my $more=shift;
	my $isAry=0;
	my @mores=();
    $isAry=1 if !defined($col);
    my $lastIdx;
	my @cols;
	if (defined($col)) {
      @cols=split(/:/,$col);
	}	

	#给定始终区间
	my $usesub=0;
	my $firstIdx=0;
	if ($more=~/:/) {
	  my @xx=split(/:/,$more);
	  $lastIdx=$xx[1];
	  $firstIdx=$xx[0];
	  $usesub=1;
	}

	if ($isAry and !$usesub) {
      $lastIdx=@$mtx-1;
	} elsif (!$isAry and !$usesub) {
	  $lastIdx=$#$mtx;
	}
	return [] if $lastIdx==-1;

	die "getIntervals(): error more" if $lastIdx<$firstIdx;

	my ($cur,$startIdx,$endIdx,$ret);
    if ($isAry) {
	  $cur=$$mtx[0];
	  $endIdx=$startIdx=$firstIdx;
	  if ($lastIdx==0) {
		push(@{$ret},[0,0]);
		push(@mores,$cur);
	  } else {
		for my $i(($firstIdx+1)..$lastIdx) {
		  if ($$mtx[$i] eq $cur) {
			$endIdx++;
		  } else {
			push(@{$ret},[$startIdx,$endIdx]);
			push(@mores,$cur);
			$startIdx=$endIdx+1;
			$endIdx++;
			$cur=$$mtx[$i];
		  }
		  if ($i==$lastIdx) {
			push(@{$ret},[$startIdx,$lastIdx]);
			push(@mores,$cur);
			last;
		  }
		} #for i
     }
   } #isAry
   else {
	$cur=join(",",@{$mtx->[0]}[@cols]);
	$endIdx=$startIdx=$firstIdx;
	if ($lastIdx==0) {
	  push(@{$ret},[$startIdx,$lastIdx]);
	  push(@mores,$cur);
	} else {
	  for my $i (($firstIdx+1)..$lastIdx) {
		if (join(",",@{$mtx->[$i]}[@cols]) eq $cur) {    
			$endIdx++;  			
		}  
		else {	  
		  push(@{$ret},[$startIdx,$endIdx]);
		  push(@mores,$cur);
		  $startIdx=$endIdx+1;
		  $endIdx++; 
		  $cur=join(",",@{$mtx->[$i]}[@cols]);	
		}
		if ($i==$lastIdx) {
		  push(@{$ret},[$startIdx,$lastIdx]);
		  push(@mores,$cur);
		  last;
		}
	  }	#for $i	  
	}
   }

   if ($more) {
	 die "getIntervals: rows in mores not same as in ret." if $#mores!=$#$ret;
	 for my $i(0..$#$ret) {
	   push(@{$ret->[$i]},$mores[$i]);
	 }
   }
  return($ret);
}

#############################################################################
#  fillGaps($mtx,$colstart,$colend):$idxmtx (prevIdx,nextIdx)
#  注意: mtx是排序后的!!!!
#  填充空隙，比如gene间填充igt
#  mtx: 已按colstart列排序的矩阵
#  colstart,colend: start/end坐标的列号（首0）
#  输出idxmtx: (prevIdx,nextIdx)表示mtx中的startIdx行和endIdx行是该GAP的首末，即GAP的长度是 mtx.nextIdx.start-mtx.prevIdx.end+1
#Usage: 
#my $ret=fillGaps($mtx,1,2)

#测试
#my $ret=fillGaps($mtx,0,1);
#my $ipre,$inext;
#if ($ret->[0][0]==-1) {
#  $inext=$ret->[0][1];
#  print "1\t".($mtx->[$inext][0]-1)."\n";
#}
#for my $i (1..($#$ret-1)) {
#  $ipre=$ret->[$i][0];
#  $inext=$ret->[$i][1];
#  print ($mtx->[$ipre][1]+1)."\t".($mtx->[$inext][0]-1)."\n";
#}
#$ipre=$ret->[$#$ret][0];
#print ($mtx->[$ipre][1]+1)."\t theEnd \n";
##############################################################################
sub fillGaps {
	my($mtx)=shift;
	my $cs=shift;
	my $ce=shift;
	my ($i,$pre,$next,$ret,$ovp);
	if ($mtx->[0][$cs]>1) {
	  push(@{$ret},[-1,0]);
	}	
	$i=0;
    while ($i<$#$mtx) {
	  $ovp=0;
	  $pre=$i;
	  $next=$pre+1;
	  while ($next<=$#$mtx) {
		if ($mtx->[$next][$cs]<=$mtx->[$pre][$ce] and $mtx->[$next][$cs]>=$mtx->[$pre][$cs]) {
			$ovp=1;
			if ($mtx->[$pre][$ce]>$mtx->[$next][$ce]) {
			  $next++;
			} else {
			  $pre=$next;
			  $next++;
			}
		} else {
		  $ovp=0;
		  last;
		}
	  }
	  #print "pre=$pre,next=$next\n";
	  $next=$#$mtx if $next>$#$mtx;
	  if (!$ovp and $next>$pre) {
		push(@{$ret},[$pre,$next]);
	  }
	  $i=$next;
    }
	if ($ovp) {
		push(@{$ret},[$pre,-1]);
	} else {
		push(@{$ret},[$next,-1]);  
	}
	return $ret;
}

#############################################################################
## 得到连续的区域，用于doAMB，做法参照IRanges.disjoin()
## $outmtx=disjoin($inmtx,$startcol=0,$endcol=1,$startrow=0,$endrow=lastRow)
## inmtx是>=两列矩阵，指定start/end列，也可以指定inmtx的哪些行，输出outmtx是按start排序后的
#############################################################################
sub disjoin {
	my ($g1,$si,$ei,$ri,$re)=@_;
	$si=0 if !defined($si);
	$ei=1 if !defined($ei);
	$ri=0 if !defined($ri);
	$re=$#$g1 if !defined($re);
	my @starts=();
	my @ends=();
	my @starts1=();
	my @ends1=();
	my $g2=[];
	for my $i($ri..$re) {
	  push(@starts,$g1->[$i][$si]);
	  push(@ends,$g1->[$i][$ei]);
	  push(@starts1,$g1->[$i][$si]-1);
	  push(@ends1,$g1->[$i][$ei]+1);
	  push(@$g2,[($g1->[$i][$si],$g1->[$i][$ei])]);
	}

    #print "\n".mtx2str($g2)."\n\n";

	@starts=uniq @starts;
	@ends=uniq @ends;

	push(@starts,@ends1);
	@starts=uniq @starts;
	#my @adjstart=sort(@starts); !!错
	my @adjstart = sort { $a <=> $b } @starts;

	pop(@adjstart);
	push(@ends,@starts1);
	@ends=uniq @ends;
	my @adjend = sort { $a <=> $b } @ends;
	@adjend=@adjend[1..$#adjend];
	my $adj=[];
	for my $i(0..$#adjstart) {
	  push(@$adj,[$adjstart[$i],$adjend[$i]]);
	}

	$g2 =sortInt($g2,0);
	$adj =sortInt($adj,0);

	my $idx=findOverlaps($adj,$g2,0,1,0,1,0,0,'first',1,'any',1,0);
	my $out=[];
	for my $i(0..$#$idx) {
	  push(@$out,[@{$adj->[$idx->[$i][0]]}]);
	}
	return($out);
}


#****************************************************************************
# Sequencing
#****************************************************************************

#############################################################################
#  seqFormat($seqfile):fa/fq/''
#  useage: $format=seqFormat($file)
#  说明: 判断序列文件是fa或fq格式,只通过第1行的标题是>还是@来判断格式
#############################################################################
sub seqFormat {
  my $f=shift;
  open(XX,"<$f") or return '';
  my $line=<XX>;
  return('fa') if (substr($line,0,1) eq '>');
  return('fq') if (substr($line,0,1) eq '@');
  return('');
}

#############################################################################
#  isBad($seq,0.9):0=good/1=bad
#  useage: if (isBad('aaaa')).. if (isBad('aaa',0.9))
#  说明: qc,判断序列是否含 10%的N 或 QT%的ATCG..
#############################################################################
sub isBad {
 my $aseq=shift;
 my $QT=0.8;
 $QT=shift if @_>0;
 my $ls=length($aseq);
 return 1 if $ls==0;
 my $an = $aseq =~ tr/Aa/Aa/;
 my $tn = $aseq =~ tr/TtUu/TtUu/;
 my $cn = $aseq =~ tr/Cc/Cc/;
 my $gn = $aseq =~ tr/Gg/Gg/;
 my $nn = 0;
 ++$nn while($aseq =~/[^ATCGU]/gi);
 if ($an/$ls>=$QT or $tn/$ls>=$QT or $cn/$ls>=$QT or $gn/$ls>=$QT or $nn/$ls>=0.1) {
   return(1);
 }
 return(0);
}

#############################################################################
## remove12($file,$sn); remove12($file,$sn,$of); 
## 去除sn列指定的seq_name列的/1/2tag,重写原文件
## sn: seq_name所在列,首1. of: 输出,若未指定,则重写原文件
## 会自动根据第1行判断是否有1/2tag,若有则重写,无则直接返回
#############################################################################
sub remove12 {
  my ($file,$sn,$of)=@_;   
  my ($line);
  $sn--;
  die "sn at least 1" if $sn<0;
  ##自动根据第1行判断是否为/1/2
  open(ININ,"<$file") or die "Cannot read $file";
  $line=<ININ>;
  $line=trim($line);
  my @nn=split(/\t/,$line);
  my $last2=substr($nn[$sn],length($nn[$sn])-2,2);
  if ($last2 ne '/1' and $last2 ne '/2') {
	close(ININ);
	return(0);
  }
  close(ININ);

  open(ININ,"<$file") or die "Cannot read $file";
  my $f2;
  if (!$of or $of eq $file) {
	$f2="$file.XXXremove12";
  } else {
    $f2=$of;
  }
  open(OO,">$f2") or die "Cannot write $f2";
  while ($line=<ININ>) {
	$line=trim($line);
	next if $line eq '';
	@nn=split(/\t/,$line);
    $nn[$sn]=substr($nn[$sn],0,length($nn[$sn])-2);
	print OO join("\t",@nn)."\n";
  }
  close(OO);
  close(ININ);
  rename "$f2","$file" if (!$of or $of eq $file);
  return(1);
}

#############################################################################
## @files=splitRaw($file,$max,$cnt,$suf,$remove12,$seqfld,$idx);
## 根据seq_name列,将file划分成小文件,返回划分后的文件列表
## max/cnt: 用于确定分组的最大值及组数
## suf: 文件后缀如.part,则输出.part1..N+1
## remove12: 是否删除seq_name后的/1/2(自动判断有无)
## seqfld: seq_name所在的列,首1. 默认1,第1列.
## seqfld类似: MCIC-SOLEXA:2:X:529:1316#0/1,根据X拆分
## idx:  首0,默认为2（适于MCIC的情况），指定seqname可以识别用于split的部分
## ??@HWI-ST741:189:C0GU5ACXX:8:2315[这列是变的]:21296:100820 1:N:0: 这种咋办
## 输出$file.part1..6文件名,各组表示从1~20,21~40..81~100,>100
## Ex. @files=splitRaw($file,100,5,'.part',1);
## 2012-06-03 默认情况下idx=2，对于LI16G的item是4，最大可到2500左右
##  @files=splitRaw($file,2500,50,'.part',1,4);
#############################################################################
sub splitRaw {
  my ($file,$max,$grpcnt,$suf,$rm,$sn,$sidx)=@_;
  die "Error grpcnt or max" if (!$grpcnt or !$max);
  $sidx=2 if $sidx eq '';
  if (!$sn) {
	$sn=0;
  } else {
	$sn--;
  }
  my (@hs,@fs,$f);
  $rm=0 if !$rm;
  
  #如果remove/1/2,则先判断是否含,只根据第1行判断
  if ($rm) {
	open(RAW,"<$file") or die "Cannot read $file!";
	my $line=<RAW>;
	my @ll=split(/\t/,$line);
    my $last2=substr($ll[$sn],length($ll[$sn])-2,2);
    if ($last2 eq '/1' or $last2 eq '/2') {
	  $rm=1;
	} else {
	  $rm=0;
	}
	close(RAW);
  }

  open(RAW,"<$file") or die "Cannot read $file!";
  #建输出文件句柄
  for my $i(1..($grpcnt+1)) {
	my $h;
    $f="$file$suf$i";
	push(@fs,$f);
	open($h,">$f");
	push(@hs,$h);
  }
  my $last=$#hs;

  #分到不同文件
  my $cnt=int($max/$grpcnt);
  $max=$grpcnt*$cnt;
  my ($line,@cols,@item,$idx,$name,$mod,$grp);
  while ($line=<RAW>) {
	$line=trim($line);
	next if $line eq '';
	@cols=split(/\t/,$line);
	if ($rm) {
      $cols[$sn]=substr($cols[$sn],0,length($cols[$sn])-2);
	}
	#确定所有的组
	@item=split(/:/,$cols[$sn]);
	if ($#item<$sidx or $item[$sidx]>$max or $item[$sidx]<=0) {
	  my $HF=$hs[$last];
      print $HF join("\t",@cols)."\n";
	  next;
	}
	$mod=$item[$sidx]%$cnt;
	$grp=int($item[$sidx]/$cnt);
	$grp-- if $mod==0;
	#print "$name, i2=$item[$idx] grp=$grp hs=$#hs\n";
	my $HF=$hs[$grp];
	print $HF join("\t",@cols)."\n";
  }

  for my $i(0..$#hs) {
	close($hs[$i]);
  }
  return(@fs);
}



#############################################################################
## @uniqs=getUniqByCols($file,$cols,$sep='');
## 根据cols指定的列,取得uniq值，用sep连接
## sep默认为''
## Ex. @uniqs=getUniqByCols($file,'0:1:2','-');
#############################################################################
sub getUniqByCols {
  my ($file,$acols,$sep)=@_;
  $acols=0 if !$acols;
  my @cols=split(/:/,$acols);
  $sep='' if !$sep;
  
  open(GUB,"<$file") or die "getUniqByCols: Cannot read $file!";
  my ($line,@items,@vals);
  while ($line=<GUB>) {
	$line=trim($line);
	next if $line eq '';
	@items=split(/\t/,$line);
	push(@vals,join($sep,@items[@cols]));
  }
  close(GUB);
  return (uniq @vals);
}


#############################################################################
## $files=splitFileByCols($file,$cols,$lbl);
## 根据cols指定的列,将file划分成小文件,返回划分后的文件列表-file.spartX及对应的cols值
## 输出$files=[file.spart1,chr1**+] 用**隔开cols
## 主要考虑到cols可能含有非法的文件名，所以用.spartX输出文件
## Ex. $files[filename,uniq]=splitFileByCols($file,'0:1:2');
##     则输出为 file.xx1..xxN
##     $files[filename,uniq]=splitFileByCols($file,'0:1:2','xx');
## 2012-10-26 增加nochk选项，若为1,则不检验uniqs个数 (用于Mtr.）
#  splitFileByCols($pairfile,"$tanchr2:$tanstrand2",'',1); 
#############################################################################
sub splitFileByCols {
  my ($file,$acols,$lbl,$nochk)=@_;
  $lbl='spart' if !$lbl;
  $acols=0 if !$acols;
  my @cols=split(/:/,$acols);
  $nochk=0 if !$nochk;
  #得到uniq值
  my (%hs,$fs,$f);
  my @uniqs=getUniqByCols($file,$acols,'**');
  if (!$nochk) {
   die "splitFileByCols: uniqs more than 200; please check $acols!!" if $#uniqs>=200;
 }

 my $nhandle=scalar(@uniqs);
 my $max=1000; #这里用3万又会出错，而实际试过32267又可以，只能设置1000了
 my $ngrp=int($nhandle/$max);
 my $rest=$nhandle%$max;
 if ($nhandle%$max!=0) {
   $ngrp+=1;
 }

for my $ng(1..$ngrp) {
  %hs=();
  my $uqs=$max*($ng-1);
  my $uqe=$uqs+$max-1;
  if ($uqe>$nhandle-1) {
	$uqe=$nhandle-1;
  }
  #print "group($ngrp): $uqs..$uqe\n";
  for my $i($uqs..$uqe) { ##2012-10-26 当打开的文件很多（比如有32845个，32267个还是可以的）会出现文件不能打开的错误，即句柄开了太多了
	my $h;
    $f="$file.$lbl.$i";
	push(@{$fs},[$f,$uniqs[$i]]);
	open($h,">$f");
	$hs{$uniqs[$i]}=$h; #chr**+对应h句柄
  }

  open(SFC,"<$file") or die "splitFileByCols: Cannot read $file! Total ".scalar(@uniqs)." uniqs\n";
  #分到不同文件
  my ($line,@items,$val);
  while ($line=<SFC>) {
	$line=trim($line);
	next if $line eq '';
	@items=split(/\t/,$line);
	$val=join('**',@items[@cols]);
	if (defined($hs{$val})) {
	  my $HF=$hs{$val};
	  print $HF "$line\n";
	}
  }
  close(SFC);
  for my $i($uqs..$uqe) {
    close($hs{$uniqs[$i]});	
  }
} #ng
  return($fs);
}

#****************************************************************************
# APA
#****************************************************************************
#############################################################################
#  findTail($seq,$from,$AT,$gap):$palen
#  useage: $palen=findTail($seq,100,'AAAAAA',4)
#  说明: 从seq的from(不包含from位置,首1)开始连续的找A或T,允许首个A/T离from的4个位置以内. 
#  如 12345AAAAAAAAAAAA from=1,gap=4,则返回扩展后的palen,若找不到palen则返回0
#  如$from=-100,则从from的左边开始匹配. AA34567 from=-3 则从3旁边的第1个A开始匹配
#############################################################################
sub findTail {
my($seq,$from,$AT,$gap)=@_;
my($palen)=0;
my($idx,$i,$nt);
$seq=uc($seq);

if ($from>0) {
	$idx=index($seq,$AT,$from); #idx=5
	if ($idx==-1 or $idx-$from>$gap) {
	  return($palen);
	}
	$nt=substr($AT,0,1);
	$palen=length($AT);
	for $i(($idx+length($AT))..length($seq)) {
	  if (substr($seq,$i,1) eq $nt) {
		$palen++;
	  }else{
		return($palen);
	  }
	}
}else{
	$from=abs($from);
	my($s)=substr($seq,0,$from-1);
	my($len)=length($s);
	$idx=rindex($s,$AT,$len); 
	#my($ii)=$len-$idx-length($AT);
	#print "idx=$idx ii=$ii\n";
	if ($idx==-1 or $len-$idx-length($AT)>$gap) {
	  return($palen);
	}
	$nt=substr($AT,0,1);
	$palen=length($AT);
	for ($i=$idx-1;$i>=0;$i--) {
	  if (substr($s,$i,1) eq $nt) {
		$palen++;
	  }else{
		return($palen);
	  }
	}
}

return($palen);
}


#****************************************************************************
# PAT
#****************************************************************************
#############################################################################
#  trimseq($longseq,$center,$left,$right,$isRC):$str
#  useage: $str=trimseq(\$longseq,301,300,99,1);        
#  说明: 传longseq的引用,$isRC=1需要反转互补,0不需要
#############################################################################
# str,center,left,right
sub trimseq {
  my ($longseq,$center,$left,$right,$isRC)=@_;
  my ($len,$offset,$subseq);
  $offset=$center-$left-1;
  $len=$left+$right+1;
  $subseq=substr($$longseq,$offset,$len);
  return $subseq if !$isRC;
  #反转互补
  return reverseAndComplement($subseq,1,1);
}

#############################################################################
#  reverseAndComplement($seq,$r,$c):$str
#  useage: $str=reverseAndComplement($seq,1,1);        
#  说明: $r=1 reverse $c=1 complement
#############################################################################
sub reverseAndComplement {
  my ($seq,$r,$c)=@_;
  return $seq if (!$r and !$c);

  $seq=uc(trim($seq));
  if ($c and $r) {
	$seq=~tr/ATCGUatcgu/TAGCATAGCA/;
	my $tmp=reverse($seq); #尼马，这里要赋给一个值，不然调用会出错
	return $tmp;
  }

  if ($c and !$r) {
	$seq=~tr/ATCGUatcgu/TAGCATAGCA/;
	return $seq;
  }

  if ($r and !$c) {
	my $tmp=reverse($seq);
	return $tmp;
  }

}


#############################################################################
#  reverseAndComplement($seq,$r,$c):$str
#  useage: $str=reverseAndComplement($seq,1,1);        
#  说明: $r=1 reverse $c=1 complement
#############################################################################
sub reverseAndComplement2 {
  my ($seq,$r,$c)=@_;
  return $seq if (!$r and !$c);
  #反转互补
  $seq=uc(trim($seq));
  my (@s)=split(//,$seq);
  my ($s2,$i);
  if ($r and $c) {
	  for ($i=$#s; $i>=0 ; $i--) {
		if ($s[$i] eq 'A') {
		  $s2.='T';
		} elsif ($s[$i] eq 'T' or $s[$i] eq 'U') {
		  $s2.='A';
		} elsif ($s[$i] eq 'G') {
		  $s2.='C';
		} elsif ($s[$i] eq 'C') {
		  $s2.='G';
		} else {
		  $s2.='N';
		}
	  } #for
  }elsif($r) {
	  for ($i=$#s; $i>=0 ; $i--) {
		$s2.=$s[$i];
	  }
  }elsif($c) {
	  for ($i=0; $i<=$#s ; $i++) {
		if ($s[$i] eq 'A') {
		  $s2.='T';
		} elsif ($s[$i] eq 'T' or $s[$i] eq 'U') {
		  $s2.='A';
		} elsif ($s[$i] eq 'G') {
		  $s2.='C';
		} elsif ($s[$i] eq 'C') {
		  $s2.='G';
		} else {
		  $s2.='N';
		}
	  } #for
  }
  return $s2;
}

#############################################################################
#  dna2amino($dna,$shift):$amino
#  useage: $str=dna2amino($dna,1);        
#  说明: shift=1(默认),2,3; 1表示从第1个位置开始翻译
#  返回amino类似：aNDRSATISYKNPGATIIg 会包括前后没翻译的核苷(输入DNA若是大写，自动转为小写)
#############################################################################
sub dna2amino {
  my($dna,$shift)=@_;
  if (!$shift) {
	$shift=1;
  }
  $dna=lc($dna);
  my $scrap = substr($dna,0,$shift-1);
  my $main = substr($dna,$shift-1);
  $main =~ s/(...)/"$convertor{uc $1}" || "?"/eg;
  return $scrap.$main;
}

#############################################################################
#  trypsin($amino,$notJunc):$amino
#  useage: $str=trypsin($amino,0);        
#  说明: notJunc:0只返回最后的一段(默认)；1返回全部，以；隔开R和K
# 例子：
#MACTWGKPELPHEASRTGHECNVKKKKKK
#MACTWGKPELPHEASR;TGHECNVK;K;K;K;K;K; （notJunc=1时）
#0.MACTWGKPELPHEASR
#1.TGHECNVK
#2~6. K
#最后是返回 TGHECNVK （notJunc=0）
#############################################################################
sub trypsin {
my ($text,$notJunc)=@_;
$notJunc=0 if !$notJunc;

my @str=split(//,$text);
my $dot='';

for my $i(0..($#str-1)) {
  if ($str[$i]=~/K|R/ and $str[$i+1] ne 'P') {
	 $dot.=($str[$i].';');
  } else {
	 $dot.=$str[$i];
  }
}
if ($str[$#str]=~/K|R/) {
  $dot.=$str[$#str].';';
}

if ($notJunc) {
  return($dot);
}

my @items=split(/;/,$dot);
for (my $i=$#items;$i>=0;$i--) {
  if (length($items[$i])>1) {
	return($items[$i]);
  }
}
return('');

}


#############################################################################
#  trimseq($longseq,$from,$to,$isRC):$str
#  useage: $str=trimseq(\$longseq,270,300,0);        
#  说明: 传longseq的引用,$isRC=1需要反转互补,0不需要
#############################################################################
# str,center,left,right
sub trimseqFromTo {
  my ($longseq,$from,$to,$isRC)=@_;
  return trimseq($longseq,$from,0,$to-$from,$isRC); 
}

#############################################################################
#  grpSame() group same pos ($dist=0)
#  useage: @grp=grpSame(\@coord); 
#  注意传递的是数组的引用!
#############################################################################
sub grpSame {
  my ($curGrpNum,$i,$j,$N,$coord,@grp);
  $coord=$_[0];

  $N=@$coord;
  $curGrpNum=1;
  $grp[0]=1;
  for ($i=1;$i<=$N-1;$i=$i+1) {
    if ($$coord[$i]==$$coord[$i-1]) {
	  $grp[$i]=$curGrpNum;
    } else {
	  $curGrpNum++;
      $grp[$i]=$curGrpNum;
	}
  }
  return @grp;
}

#############################################################################
#  grpByPos() group by position within dist
#  useage: @grp=grpByPos(\@coord,$dist); 
#  注意传递的是数组的引用!
#  计算两两差,N个数,计算N-1组diff
#  遍历diff,直到sum(diff)>dist,则要去除头或尾的一个
#    去头:第1个设置组号,其它的继续放回等待
#    去尾:除最后1个继续放回等待,前面的设置组号
#############################################################################
sub grpByPos_I {
  my ($dist,$curGrpNum,$i,$j,$N,$sumfrom,$sum,$coord,@diff,@grp);
  $coord=$_[0];
  $dist=$_[1];

  #print "coord=@$coord","\n","dist=$dist\n";
  
  $N=@$coord;
  $curGrpNum=1;
  #只有1个元素
  if ($N==1) {
  	$grp[0]=$curGrpNum;
  
  } else { #2+个元素
    $sum=0;
    $sumfrom=0;
    for ($i=0;$i<$N-1;$i=$i+1) {
      $diff[$i]=$$coord[$i+1]-$$coord[$i];
      $sum=$sum+$diff[$i];
      if ($sum>$dist) {
		      #若此分组最后1个与上一分组最后1个位置距离<=dist,则合并到前一组
			  if ($$coord[$i]-$$coord[$sumfrom-1]<=$dist) {
				$curGrpNum--;
              }
    	      for ($j=$sumfrom;$j<=$i;$j++) {
    			$grp[$j]=$curGrpNum;
    	      }
    		  $curGrpNum+=1;
    		  $sumfrom=$i+1;
    		  if ($sumfrom==$N-1) { #只剩1个
				  if ($$coord[$sumfrom]-$$coord[$sumfrom-1]<=$dist) {
					$curGrpNum--;
				  }
    			$grp[$sumfrom]=$curGrpNum;
    			$curGrpNum++; 
				$sumfrom++;
    		  }
    		  $sum=0;
      }
    } #for
  } #else
  
  #剩下的dist内的所有元素
  #若最后1组的最后1个也与前1组最后1个距离<=dist,则合并到前一组
  if ($$coord[$N-1]-$$coord[$sumfrom-1]<=$dist) {
    $curGrpNum--;
  }
    for ($j=$sumfrom;$j<=$N-1;$j++) {
      $grp[$j]=$curGrpNum;
    }
  
  #print "diff=@diff","\n";
  #print "grp=@grp","\n";
  return @grp;
}

#############################################################################
#  grpByPos() group by position within dist
#  useage: @grp=grpByPos(\@coord,$dist); 
#  注意传递的是数组的引用!
#  比较前后2个位置，如果距离<=dist,则合并,不管合并后的宽度
#############################################################################
sub grpByPos {
  my $coord=$_[0];
  my $dist=$_[1];
  
  my $N=@$coord-1;
  return () if $N<0;
  my (@grp);
  #只有1个元素
  if ($N==0) {
  	$grp[0]=0;  
  } else { #2+个元素
    $grp[0]=0;
	for my $i(1..$N) {
	  if ($$coord[$i]-$$coord[$i-1]<=$dist) {
		$grp[$i]=$grp[$i-1];
	  } else {
		$grp[$i]=$grp[$i-1]+1;
	  }
    }
  }
  return @grp;
}

#############################################################################
#  isIP($papos,$nt,$nts,$seq):$int
#  useage: $IP=isIP($papos,$nt,$nts,\$seq);         
#  说明:不判断边界,默认位点左右各有9,10nt (-10~-1[PA],1~10)
#  返回1代表是IP,否则0
#############################################################################
sub isIP {
	my($papos,$nt,$nts,$seq)=@_;
	my($s,$e,$ee,$cnt,$subseq,@wseq,$i);
	$s=0;
	#左取10nt(含pa),右取10nt,共20nt,下标从0开始	
    if ($papos-10<0) { #2011/3/6修改，原来是无这个判断，对开头的序列判断是错的!
	  my $start=0;
      $subseq=uc(substr($$seq,$start,10+$papos));
	} else {
	  $subseq=uc(substr($$seq,$papos-10,20));
	}
    return 1 if (index($subseq,$nts)>=0); #2011/3/6修改，原来是 >0 ，对第一个开头含As是错的!
	@wseq=split(//,$subseq);
	$e=$s+9;	
	for $i ($s..$e) {
		$cnt++ if($wseq[$i] eq $nt);
	}
    #print "$s~$e cnt=$cnt\n";
	return 1 if($cnt>=7);
	$ee=$s+19;
	$e++;	
    while ($e<=$ee) {
		$cnt++ if($wseq[$e] eq $nt);
		$cnt-- if($wseq[$s] eq $nt);
		$s++;
		$e++;
		#print "$s~$e s=@wseq[$s] e=@wseq[$e] cnt=$cnt\n";
        return 1 if($cnt==7);
    }
	return 0;
}

#############################################################################
#  formatPatOutput($file);
#  useage: formatPatOutput($file);;
#  说明:#格式化patronus的输出,并替换原文件
#############################################################################
sub formatPatOutput {
my($f)=shift;
my($nl)=1;
my($line,$s,$i1,$i2,$i3,@ss,$i,$j,$mtx,$idx);
  open(IN, "<$f")|| (print "cannot open $f\n" and return);
  $idx=0; 
  while($s=<IN>) {
	if ($s=~/^\#/) {
	  $s=~s/[^a-zA-ZuUnN]//g;
	  $mtx->[$idx][0]=$s;
	} elsif ($s=~/^P/) {
	  $s=trim($s);
	  $i1=index($s,'(');
	  $i2=index($s,')');
	  $i3=index($s,'=');
      $mtx->[$idx][1]=substr($s,$i1+1,$i2-$i1-1);
      $mtx->[$idx][2]=substr($s,$i3+1,length($s)-$i3-1);
	} elsif ($s=~/^Z/) { #Z_value(213)=4.91956 Mean=146.156 StdDev=13.5875
	  @ss=split(/\s+/,$s);
	  for $i(0..0) { #只要Z_value
		$i3=index($ss[$i],'=');
		$s=trim($s);
		$s=substr($ss[$i],$i3+1,length($ss[$i])-$i3-1);
		$mtx->[$idx][3+$i]=$s;
	  }
	} 	
	if ($nl%3==0) {
	  $idx++;
	}
   $nl++;
  }#while
  close IN;
 
 #Occu_Rank
 #moitf occurence p_value Z_score Occu_Rank PZ_Rank
 my(@occ,@order);
 for $i(0..$#$mtx) {
	push(@occ,$mtx->[$i][1]);
 }
 @order=shellSort(0,$#occ,'desc',@occ);
 for $i(0..$#order) {
	for $j(0..$#order) {
		if ($order[$j]==$i) {
		  push(@{$mtx->[$i]},$j+1);
		}
	}	
 }

#PZ_Rank
@occ=();
 for $i(0..$#$mtx) {
	push(@occ,$mtx->[$i][2]);
 }
 @order=shellSort(0,$#occ,'asc',@occ);
 for $i(0..$#order) {
	for $j(0..$#order) {
		if ($order[$j]==$i) {
		  push(@{$mtx->[$i]},$j+1);
		}
	}	
 }
my($mtx2);
push(@$mtx2,[('moitf', 'occurence', 'p_value', 'Z_score', 'Occu_Rank','PZ_Rank')]);
for $i(0..$#order) {
  push(@$mtx2,[@{$mtx->[$order[$i]]}]);
}

 #保存到文件
 saveMtx2File($mtx2,"${f}x");
 rename("${f}x", $f);
}

#############################################################################
# getNt(@prb):chr
#调用: $c=getNt(@prb);  ATCG
#说明: 根据背景概率(整数的),产生联子
#############################################################################
sub getNt{
  my (@prb)=@_;
  my ($nt,$sum,$i,$r);
  $nt='A';
  $sum=0;
  for $i (0..3) {
	$sum+=$prb[$i];
  }
  $r=int(rand($sum))+1;
  #print "sum=$sum,@prb";
  for $i (1..3) {
	$prb[$i]=$prb[$i]+$prb[$i-1];
  }
  if ($r>=1 and $r<=$prb[0]) {
	return 'A';
  } elsif ($r>=1+$prb[0] and $r<=$prb[1]) {
	return 'T';
  } elsif ($r>=1+$prb[1] and $r<=$prb[2]) {
	return 'C';
  } elsif ($r>=1+$prb[2] and $r<=$prb[3]) {
	return 'G';
  } else {
	die "getNt() No ATCG!\n";
  }
}


#****************************************************************************
# Sequences / PAS
#****************************************************************************
#############################################################################
#  countSubstr($str,$substr):$count
#  useage: $i=countSubstr('AAATTXXXX','AT');      
#  说明: 计算字符串中给定字符的个数
#############################################################################
sub countSubstr {
my ($str,$substr)=@_;
my $cc = 0;
my $tmp = 0;
if( $tmp = () = ($str =~ /$substr/g ) ) { 
	$cc += $tmp;
}
return $cc;
}

#############################################################################
#  convertChrToIdx($chr):$idx
#  useage: $i=convertChrToIdx('A');      
#  说明: atcgATCG-->下标
#############################################################################
sub convertChrToIdx {
  my($i)=uc(shift);
  if ($i eq 'A') {
	return 0;
  } elsif ($i eq 'T' or $i eq 'U') {
	return 1;
  } elsif ($i eq 'C') {
	return 2;
  } elsif ($i eq 'G') {
	return 3;
  } 
}

#############################################################################
#  convertIdxToChr($idx):$chr
#  useage: $c=convertIdxToChr(0);      
#  说明:下标-->ATCG
#############################################################################
sub convertIdxToChr {
  my($i)=shift;
  if ($i==0) {
	return 'A';
  } elsif ($i==1) {
	return 'T';
  } elsif ($i==2) {
	return 'C';
  } elsif ($i==3) {
	return 'G';
  } 
}

#############################################################################
#  getKgramId($gram):$i
#  useage: $id=getKgramId('AATAAA');      
#  说明:根据kgram,得到下标
#############################################################################
sub getKgramId {
  my $g=shift;
  if ($g=~/N/) {
	return(-1); #不判断带N的联子
  }
  my(@gram)=split(//,$g);
  my($k)=$#gram+1;
  my ($i,$idx,$rt);
  $rt=0;
  for $i(1..$k){
    $idx=convertChrToIdx($gram[$i-1]);
	next if $idx==0;
	$rt+=$idx*(4**($k-$i));
  }
  return $rt;
  exit 0;
}

#############################################################################
#  genOneKGram($k,$idx):$gram
#  useage: $gram=genOneKGram(6,0);      
#  说明:产生1个空的k联子，idx表示下标
#############################################################################
sub genOneKGram {
  my($k)=shift;
  my($idx)=shift;
  my ($s,$j,$res,$num);
    $num=$idx;
	for $j(1..$k) {
	  $res=$num%4;
	  $num=int($num/4);
	  $s=convertIdxToChr($res).$s;
	}
  return $s;
}

#############################################################################
#  genKgrams($k,$withvalue):@kgrams
#  useage: @kgrams=genKgrams(6,0);      
#  说明:产生空的k联子，withvalue表示是否右边产生＝0
#############################################################################
sub genKgrams {
  my($k)=shift;
  my($v)=shift;
  my ($cnt,$i,$num,$s,$j,$res,@kgram);
  $s='';
  $cnt=4**$k;
  for $i(0..$cnt-1){
    $num=$i;
	for $j(1..$k) {
	  $res=$num%4;
	  $num=int($num/4);
	  $s=convertIdxToChr($res).$s;
	}
	$s.=('=0') if $v==1;
    push(@kgram,$s);
	$s='';
  }
  return @kgram;
}

#############################################################################
#  cntKgramsByK($seqfile,$from,$to,$k,$gapOronce):@cnts
#  useage: 
#  常规滑动方式($gapOronce不设置或=0): @cnts=cntKgrams('xx.fasta',-1,-1,6); 
#  gap方式($gapOronce>0): @cnts=cntKgrams('xx.fasta',-1,-1,6,2); 
#  once方式($gapOronce<0): @cnts=cntKgrams('xx.fasta',-1,-1,6,-1);      
#  说明:统计全部k-gram在seqfile中出现次数,$from,$to<1则为统计整条序列
#  2016/2/21 add w=1/0 weighted by seqfile title PAT number
#############################################################################
sub cntKgramsByK {
  my($seqfile,$from,$to,$k,$go,$w)=@_;
  my($i,$j,@cnts,$line,$s,$e,$idx,%ks);
  my($all)=0;
  $go=0 if !$go;
  $all=1 if($from<1 or $to<1);
  $s=$from-1;
  $e=$to-$k;
  open (INPUT,"<$seqfile") or die "cannot open file $_!\n";
  for $i(0..4**$k-1) {
	$cnts[$i]=0;
  }

 my $PAT=1;
 $w=0 if !defined($w);

  my $e1=$e;
  while($line=<INPUT>) {
	$line=trim($line);
	if($line=~/^>/) {
	  if ($w) {
	    my @x=split(';',$line);
		$PAT=$x[$#x];
	  }
	  next;
	}

    $e1=$e;
	if ($all==1) {
	  $s=0;
	  $e1=length($line)-$k;
	} else { #若序列太短,则取全长
	  $e1=length($line)-$k if $to>length($line);
	}

	if ($go==0) { #常规方式
	  for $i($s..$e1) {
	    $idx=getKgramId(substr($line,$i,$k));
	    $cnts[$idx]+=$PAT if $idx!=-1;
	  }
	}
	
	elsif ($go<0) { #once方式
	  %ks=();
	  $i=$s;
	  while ($i<=$e1) {
	    $idx=getKgramId(substr($line,$i,$k));
	    $ks{$idx}=$i-$s if $idx!=-1; #记录kgram和位置
	    $i++;
	  }
	  #累加到总体中
	  foreach $idx (keys(%ks)){
        $cnts[$idx]+=$PAT;
      }	    
    }
	
	else{ #gap方式
	  $i=$s;
	  while ($i<=$e1) {
	    #若gap=N,则一次取N+1组kgram判断
	    %ks=();
	    for $j(0..$go) {
		 if ($i<=$e1) {
           $idx=getKgramId(substr($line,$i,$k));		 
		   $ks{$idx}=$i-$s if $idx!=-1; #gap内只计最后位置,不增加计数	
	       $i++;	
		 }
	    }
	    foreach $idx (keys(%ks)){
          $cnts[$idx]+=$PAT  if $idx!=-1;
        }
	  }	 
	} #else

  } #while line
  close INPUT;
  return @cnts;
}


#############################################################################
#  cntKgrams($seqfile,$from,$to,$gapOronce,@grams):@cnts
#  useage: @cnts=cntKgrams('xx.fasta',-1,-1,('aataaa','atcg'));      
#  说明:统计给定grams在seqfile中出现次数,$from,$to<1则为统计整条序列
#  2016/2/21 加grams判断，末尾为1，则表示加权
#############################################################################
sub cntKgrams {
  my($seqfile,$from,$to,$go,@grams)=@_;
  my($i,$j,$k,@cnts,$line,$s,$e,$klen);
  my $all=0;
  ($all)=1 if($from<1 or $to<1);
  $go=0 if !defined($go);
  $s=$from-1;
  open (INPUT,"<$seqfile") or die "cannot open file $_!\n";
  for $i(0..$#grams) {
	$cnts[$i]=0;
  }

  my $w=0;
  if ($grams[$#grams]==1) {
	pop(@grams);
	$w=1;
  }
  my $PAT=1;

  while($line=<INPUT>){	
	$line=trim($line);
	if($line=~/^>/) {
	  if ($w) {
	    my @x=split(';',$line);
		$PAT=$x[$#x];
	  }
	  next;
	}

	#对于每个gram,单独计算(因为允许gram长度不一)
	for $k(0..$#grams) {
		#滑动窗口首尾
		$klen=length($grams[$k]);
		$e=$to-$klen;
		if ($all==1) {
		  $s=0;
		  $e=length($line)-$klen;
		} else { #若序列太短,则取全长
		  $e=length($line)-$klen if $to>length($line);
		}

		if ($go==0) { #常规方式
		  for $i($s..$e) {
			$cnts[$k]+=$PAT if substr($line,$i,$klen) eq $grams[$k];
		  }
		}
		
		elsif ($go<0) { #once方式
		  for $i($s..$e) {
			if (substr($line,$i,$klen) eq $grams[$k]) {
			  $cnts[$k]+=$PAT ;
			  last; #只计一次
			}
		  }   
		}
		
		else{ #gap方式
		  $i=$s;
		  while ($i<=$e) {
			if (substr($line,$i,$klen) eq $grams[$k]) {
			  $cnts[$k]+=$PAT ;
			  $i+=($go+1); #跳过$gap+1
			} else {
              $i++;
		    }
		  }	 
		} #else

	} #k
  } #while
  close INPUT;
  return @cnts;
}

#############################################################################
#  cntEachPosByK($seqfile,$from,$to,$k,$go):$cnts
#  useage: 
#    1) $cnts=cntEachPosByGap('xx.fasta',1,400,6); #常规方式 go=0   
#    2) $cnts=cntEachPosByGap('xx.fasta',1,400,6,3); #gap方式 go>0
#    3) $cnts=cntEachPosByGap('xx.fasta',1,400,6,-1); #once方式 go<0
#  说明:
#    1) 统计seqfile中从from~to的位置的各kgram出现次数,必须设置from,to
#    2) 若go>0,则下一个**相同**kgram与前一个相差gap个位置
#    3) 若go<0(once),则同一种kgram在同一条序列中只计一次,位置取最后出现的那个位置
#  2016/2/21 增加w选项，若1则用seqfile标题中；最后一字段作为PAT权重
#############################################################################
sub cntEachPosByK {
  my($seqfile,$from,$to,$k,$go,$w)=@_;
  my($i,$j,$cnts,$line,$s,$e,$idx,$oidx,%onceks,%gapks);
  $go=0 if !$go;
  die "From or To wrong!" if($from<1 or $to<1);
  $s=$from-1;
  $e=$to-$k;
  open (K_INPUT,"<$seqfile") or die "cannot open file $_!\n";
  #初始化个数阵
  for $i(0..4**$k-1) {
	for $j($s..$e) {
	  $cnts->[$i][$j-$s]=0;
	}	
  }
 
 my $PAT=1;
 $w=0 if !defined($w);
 my $e1;
  #计数1:常规方式
 if ($go==0) {
  while($line=<K_INPUT>){	
	$line=trim($line);
	if($line=~/^>/) {
	  if ($w) {
	    my @x=split(';',$line);
		$PAT=$x[$#x];
	  }
	  next;
	}

	#若序列太短,则取全长
	$e1=$e; #2011/12/1 BUG修正
	$e1=length($line)-$k if $to>length($line);
	$i=$s;
	while ($i<=$e1) {
	  $idx=getKgramId(substr($line,$i,$k));
	  $cnts->[$idx][$i-$s]+=$PAT  if $idx!=-1;
	  $i++;
	}
  } 
}
  #计数2:once方式:1个kgram在1条序列仅计数一次
 elsif ($go<0) {
  while($line=<K_INPUT>){	  
	$line=trim($line);
	if($line=~/^>/) {
	  if ($w) {
	    my @x=split(';',$line);
		$PAT=$x[$#x];
	  }
	  next;
	}
	%onceks=();
	#若序列太短,则取全长
	$e1=$e; #2011/12/1 BUG修正
	$e1=length($line)-$k if $to>length($line);
	$i=$s;
	while ($i<=$e1) {
	  $idx=getKgramId(substr($line,$i,$k));
	  $onceks{$idx}=$i-$s  if $idx!=-1; #记录kgram和位置
	  $i++;
	}
	#累加到总体中
	foreach $idx (keys(%onceks)){
        $cnts->[$idx][$onceks{$idx}]+=$PAT  if $idx!=-1;
    }	
  } 
 } 

 #计数3: gap方式
 elsif ($go>0) {
  while($line=<K_INPUT>){	
	$line=trim($line);
	if($line=~/^>/) {
	  if ($w) {
	    my @x=split(';',$line);
		$PAT=$x[$#x];
	  }
	  next;
	}
	#若序列太短,则取全长
	$e1=$e; #2011/12/1 BUG修正
	$e1=length($line)-$k if $to>length($line);
	#$e=$to-$k;
	$i=$s;
	#print "i=$i s=$s e=$e\t";
	while ($i<=$e1) {
	  #若gap=N,则一次取N+1组kgram判断
	  %gapks=();
	  for $j(0..$go) {
		 if ($i<=$e) {
           $idx=getKgramId(substr($line,$i,$k));		 
		   $gapks{$idx}=$i-$s  if $idx!=-1; #gap内只计最后位置,不增加计数	
		   #print substr($line,$i,$k)." idx=".$idx." pos=".$gapks{$idx}."\n";
	       $i++;	
		 }
	  }
	  foreach $idx (keys(%gapks)){
        $cnts->[$idx][$gapks{$idx}]+=$PAT if $idx!=-1;
      }
	}	
  } 
 }

  close K_INPUT;
  return $cnts;
}

#############################################################################
#  cntEachPosByGrams($seqfile,$from,$to,$go,@grams):$cnts
#  useage: $cnts=cntEachPosByGrams('xx.fasta',1,400,@grams));      
#  说明:统计seqfile中从from~to的位置的各kgram出现次数,必须设置from,to
#  2016/2/21 添加 grams中最后1行如果是1则是表示加权
#############################################################################
sub cntEachPosByGrams {
  my($seqfile,$from,$to,$go,@grams)=@_;
  my($i,$j,$cnts,$line,$s,$e,$idx,$l,%gpos,$key);
  die "From or To wrong!" if($from<1 or $to<1);
  #只能同长
  my($k)=length($grams[0]);
  $s=$from-1;
  $e=$to-$k;
  open (INPUT,"<$seqfile") or die "cannot open file $_!\n";
  for $i(0..$#grams) {
	for $j($s..$e) {
	  $cnts->[$i][$j-$s]=0;
	}	
  }
  my $w=0;
  if ($grams[$#grams]==1) {
	pop(@grams);
	$w=1;
  }
  $go=0 if !defined($go);

  my $PAT=1;
  my $e1=$e;
  while($line=<INPUT>)
	{
		$line=trim($line);
	if($line=~/^>/) {
	  if ($w) {
	    my @x=split(';',$line);
		$PAT=$x[$#x];
	  }
	  next;
	}
        #若序列太短,则取全长
	    $e1=$e; #2011/12/1 BUG修正
	    $e1=length($line)-$k if $to>length($line);
		#逐个比较
	  if ($go==0) {		 
		for $i($s..$e1) {
		  $l=substr($line,$i,$k);
		  for $j(0..$#grams) {
			 if ($l eq $grams[$j]) {
				$cnts->[$j][$i-$s]+=$PAT;
				last;
			 }            
		  }
	    }
	  }	
		elsif ($go<0){ #once 
			 %gpos=();
			 for $i($s..$e1) {
			  $l=substr($line,$i,$k);
			  for $j(0..$#grams) {
				 if ($l eq $grams[$j]) {
					$gpos{$j}=$i-$s; #用%{gramidx}=pos记录,只保存最后出现的位置			
					last;
				 }            
			  }
			 } 
			 foreach $key (keys(%gpos)) {
			   $cnts->[$key][$gpos{$key}]+=$PAT;
			 }   
		}

		elsif ($go>0){ #gap
		my($gi);
		$i=$s;
		while ($i<=$e1) {
		  #若gap=N,则一次取N+1个序列中的kgram作判断
		  %gpos=();
		  for $j(0..$go) {
			 if ($i<=$e1) {
			   $l=substr($line,$i,$k);	
			   for $gi(0..$#grams) {
				 if ($l eq $grams[$gi]) {
					$gpos{$gi}=$i-$s; #用%{gramidx}=pos记录,只保存最后出现的位置			
					last;
				 }            
			   }			   
			   $i++;	
			 }
		  } #j
		  foreach $idx (keys(%gpos)){
			$cnts->[$idx][$gpos{$idx}]+=$PAT  if $idx!=-1;
		  }
		}	
	  }

	} #while
  close INPUT;
  return $cnts;
}

#############################################################################
#  kcnt2pssm($cntmtx):$pssm
#  useage: $pssm=kcnt2pssm($kcntmtx);      
#  说明:kcnt阵转换为pssm阵
#  2010/3/5: 加入psedo_count=1
#############################################################################
sub kcnt2pssm {
  my($cntmtx)=@_;
  my($pssm,@rowsum,@colsum,$sumall,$i,$j,$p,$pg,$pk);
  #计算第N行的所有列的和
  for $i(0..$#$cntmtx) {
    for $j(0..$#{$cntmtx->[$i]}) {
      $rowsum[$i]+=$cntmtx->[$i][$j];
	}
  }
  #计算第N列的所有行的和
  for $i(0..$#{$cntmtx->[0]}) {
    for $j(0..$#$cntmtx) {
      $colsum[$i]+=$cntmtx->[$j][$i];
	}
  }
  #所有联子总和（行和或列和的和）
  for $i(0..$#colsum) {
    $sumall+=$colsum[$i];
  }

  for $i(0..$#$cntmtx) {
    for $j(0..$#{$cntmtx->[$i]}) {
     # if ($colsum[$j]==0) {
	#	$pssm->[$i][$j]=NONPSSM();
	#	next;
     # }else {
	#	$pk=$cntmtx->[$i][$j]/$colsum[$j];
	 # }
	  $pk=($cntmtx->[$i][$j]+0.25)/($colsum[$j]+1);
	  $pg=$rowsum[$i]/$sumall;
	  if ($pg==0) {
		$pg=0.00000000001;
	  }
	 # if ($pg==0) {
	#	$pssm->[$i][$j]=NONPSSM();
	#	next;
	 # }else {
     #   $p=$pk/$pg;
	 # }
	 # if ($p==0) {
	#	$pssm->[$i][$j]=NONPSSM();
	#	next;
	 # }else{
     #   $pssm->[$i][$j]=log($p)/log(2);
	 # }
	  $pssm->[$i][$j]=log($pk/$pg)/log(2);
	}
  }#for i
  return $pssm;
}

#############################################################################
#  sortKcnt($cntmtx):@order
#  useage: @order=sortKcnt($kcntmtx);      
#  说明:降序排序kcnt矩阵,以所有列的和排序
#  
#############################################################################
sub sortKcnt {
  my($cntmtx)=shift;
  my(@rowsum,$i,$j);
  #计算第N行的所有列的和
  for $i(0..$#$cntmtx) {
    for $j(0..$#{$cntmtx->[$i]}) {
      $rowsum[$i]+=$cntmtx->[$i][$j];
	}
  }
  return shellSort(0,$#rowsum,'desc',@rowsum);
}

#############################################################################
#  sortPssm($pssm):@order
#  useage: @order=sortPssm($pssm);      
#  说明:降序排序pssm矩阵,以所有列的最大值排序
#############################################################################
sub sortPssm {
  my($pssm)=shift;
  my(@rowmax,$i,$j);
  #计算第N行的所有列的max
  for $i(0..$#$pssm) {
    $rowmax[$i]=-99999999999;
  }
  for $i(0..$#$pssm) {
    for $j(0..$#{$pssm->[$i]}) {
      $rowmax[$i]=$pssm->[$i][$j] if $rowmax[$i]<$pssm->[$i][$j];
	}
  }
  return shellSort(0,$#rowmax,'desc',@rowmax);
}

#############################################################################
#  fas2tbl($fasfile):$tbl
#  useage: $tbl=fas2tbl($fasfile);      
#  说明: 读入fas,输出2列matrix(title,seq)
#############################################################################
sub fas2tbl {
  my($f)=shift;
  my($tbl);
  open (fas2tbl_INPUT,"<$f") ||  return $tbl;
  
  my($title)="";
  my($seq)="";
  my($line);
  while($line=<fas2tbl_INPUT>)
  {
	$line=trim($line);
	if($line=~/^>/)
	{
		if(length($title)>0 and length($seq)>0) {
		  push(@$tbl,[$title,$seq]);
		}
		$title=substr($line,1,length($line)-1);
		$title=~s/\t/ /g; #把TAB替换成空格,免得2列不对.
		$seq="";
	}
	else
	{
		$seq=$seq.$line;
	}
  }

if(length($title)>0 and length($seq)>0) {
  push(@$tbl,[$title,$seq]);
}
close(fas2tbl_INPUT);
return $tbl;
}


#############################################################################
#  ($poly,$polypos,$polylen,$polyedit)=findPolyAT($seq,$search=A/T/AT,tlen,tregion)
#  useage: 
#  ($poly,$polypos,$polylen,$polyedit)=findPolyAT($str,'A',8,20);
#  说明: 查找A/T，或自动判断polyA/T
#输入：
#seq=序列
#search=('A','T','AT') 若为AT，则自动查找，然后判断A或T
#tlen=8 polyA尾巴的长度（含非A部分）
#tregion=20 对polyA查找末尾20nt，对polyT查找开头20nt
#返回 $poly=A/T/N,$polypos首1,$polylen,$polyedit （若为N，则后面三个数为-1）
#例子：($poly,$polypos,$polylen,$polyedit)=findPolyAT($str,'A',8,20);

#polyA的判断：满足Aper>=0.8; Acnt>=8; Alen>=tlen 以及polyA的最后一个A要在seq尾部的tregion内
#polyT的判断：满足Tper/Tcnt，以及polyT的第一个T要在seq头部的tregion内
#正则表达式	 /findPolyAT/ 表示以8A开头，且中间可以间隔1个其它字符，但必须有AA相间隔。
#比如 AAAAAXAAXAAXAA 表示一串 而 AAAAAXAXA则不满足条件
#############################################################################
sub findPolyAT_regError { #正则表达式还是有点问题 比如 $str='AAAAAAAAXAXXXA' 用 /(A{8,}([B-Z]{0,1}AA+[B-Z]{0,1})+A+)/g 就匹配不到，或者甚至用A{6,}也不行！
my $seq=shift;
my $search=shift;
my $tlen=shift;
my $tregion=shift;

my ($polyTseq,$polyAseq,$Acnt,$Tcnt,$Aper,$Tper)=('','',0,0,0,0);
my ($polyApos,$polyAlen,$polyAedit)=(0,0,0);
my ($polyTpos,$polyTlen,$polyTedit)=(0,0,0);

my @matches=();

if ($search eq 'A' or $search eq 'AT') {
while ($seq =~ /(A{8,}([B-Z]{0,1}AA+[B-Z]{0,1})+A+)/g) { 
  push(@matches,$&); #匹配的串
  pos($seq) = pos($seq); 
} 
if ($#matches!=-1) {
  $polyAseq=$matches[$#matches]; #polyA找末尾，polyT找第一个
  $Acnt=countSubstr($polyAseq,'A');  
  $Aper=$Acnt/length($polyAseq);
  my $idx=index($seq,$polyAseq);
  $polyApos=$idx+1;
  $polyAlen=length($polyAseq);
  $polyAedit=$polyAlen-$Acnt;
  #print "($polyApos,$polyAlen,$polyAedit)\n"; 
}
}

if ($search eq 'T' or $search eq 'AT') {
@matches=();
if ($seq =~ /(T{8,}([A-SU-Z]{0,1}TT+[A-SU-Z]{0,1})+T+)/g) {  #polyT找开头，所以用if就可以
  push(@matches,$&); #匹配的串
  pos($seq) = pos($seq); 
} 
if ($#matches!=-1) {
  $polyTseq=$matches[0];
  $Tcnt=countSubstr($polyTseq,'T');  
  $Tper=$Tcnt/length($polyTseq);
  my $idx=index($seq,$polyTseq);
  $polyTpos=$idx+1;
  $polyTlen=length($polyTseq);
  $polyTedit=$polyTlen-$Tcnt;
}
}

my ($haveA,$haveT)=(($polyAseq and $polyAlen>=$tlen and $Aper>=0.8 and $Acnt>=8 and $polyApos+$polyAlen-1>=length($seq)-$tregion+1),($polyTseq and $polyTlen>=$tlen and $Tper>=0.8 and $Tcnt>=8 and $polyTpos<=$tregion));

#return ($poly,$polypos,$polylen,$polyedit);
if ( ($search eq 'A' and !$haveA) or ($search eq 'T' and !$haveT) or ($search eq 'AT' and !$haveA and !$haveT ) ) { 
  return ('N',-1,-1,-1);
} elsif ( ($search eq 'A' and $haveA) or ($search eq 'AT' and $haveA and !$haveT) or ($search eq 'AT' and $haveA and $haveT and length($polyAseq)>=length($polyTseq)) ) {
  return ('A',$polyApos,$polyAlen,$polyAedit);
} elsif ( ($search eq 'T' and $haveT) or ($search eq 'AT' and $haveT and !$haveA) or ($search eq 'AT' and $haveA and $haveT and length($polyTseq)>=length($polyAseq)) ) {
  return ('T',$polyTpos,$polyTlen,$polyTedit);
}

}


#############################################################################
#  mutateMotif(fromIdx,mutatenum,\@source,\@items,\@results)) 
#  useage: 
#  @source=split(//,'AATAAA'); @items=split(//,'ATCG'); $bits=1;
#  mutateMotif(0,$bits,\@source,\@items,\@mutes); 
#  print join("\n",@mutes);
#  说明: 最后结果在@mutes中,每个值1个motif
#############################################################################
sub mutateMotif {
  my($currIdx,$remain,$source,$items,$mutes)=@_;
  if ($remain==0) {
	push(@$mutes,join('',@$source));
	return;
  }
  if ($currIdx==scalar(@$source)) {
	return;
  }
  my $currValue=$$source[$currIdx];
  foreach my $item (@$items) {
	if ($item ne $currValue) {
	  $$source[$currIdx]=$item;
	  mutateMotif($currIdx+1,$remain-1,$source,$items,$mutes);
	}
  }
  $$source[$currIdx]=$currValue;
  mutateMotif($currIdx+1,$remain,$source,$items,$mutes);
}

#****************************************************************************
# DB / database
#****************************************************************************
#############################################################################
#  ($set,$rv)=execSql($dbh,$sql)
#  useage: ($set,$rv)=execSql($dbh,'select * from xxtbl')
#  说明:执行sql语句
#############################################################################
sub execSql {
  my($dbh,$sql)=@_;
  #print "$sql\n";
  my $sth = $dbh->prepare($sql) or die $dbh->errstr;
  my $rv=$sth->execute or die;
  my $set = $sth->fetchall_arrayref or die ; 
  $sth->finish();
  $rv=$#$set+1;
  return($set,$rv);
}

#############################################################################
#  getTableIndex($dbh,$tbl)
#  useage: $mtx=getTableIndex($dbh,'atbl')
#  说明: 取得指定表的index，返回矩阵，2列，第1列为索引名，第2列为,隔开的列，比如：index_xx chr,strand,coord
#        若表无索引，则返回空矩阵
#############################################################################
sub getTableIndex { 
  my ($dbh,$tbl)=@_;
  my $mtx=[];
	my $sql="show index from $tbl";
	my ($index,$rv)=execSql($dbh,$sql);
	return $mtx if $#$index==-1;
	my ($i,$key,$fld);
	$key=$index->[0][2];
	$fld=$index->[0][4];
	for $i ( 1 .. $#$index) { 
	  if ($key eq $index->[$i][2]) {
		$fld.=",$index->[$i][4]";
	  } else {
		push(@{$mtx},[$key,$fld]);
		$fld=$index->[$i][4];
		$key=$index->[$i][2];
	  }
	} #for $i
	push(@{$mtx},[$key,$fld]);
	return($mtx);
}



#############################################################################
#  cloneTblIndex($dbh,$fromtbl,$totbl)
#  useage: cloneTblIndex($dbh,'fromtbl','totbl')
#  说明: 在totbl中创建fromtbl的index
#############################################################################
#Table	Non_unique	Key_name	Seq_in_index	Column_name	Collation	Cardinality	Sub_part	Packed	Null	Index_type
#t_gff_org_v7	1	idx_chr	1	chr	A	9	NULL	NULL	YES	BTREE
#t_gff_org_v7	1	idx_css	1	chr	A	196	NULL	NULL	YES	BTREE
#t_gff_org_v7	1	idx_css	2	strand	A	375	NULL	NULL	YES	BTREE
#t_gff_org_v7	1	idx_css	3	start	A	347812	NULL	NULL	YES	BTREE

sub cloneTblIndex {
	my($dbh,$ftbl,$ttbl)=@_;
	my $sql="show index from $ftbl";
	my $sth = $dbh->prepare($sql) or die $dbh->errstr;
	$sth->execute or return;
	my $index = $sth->fetchall_arrayref or die ; 
	return if $#$index==-1;
	my ($i,$key,$fld);
	$key=$index->[0][2];
	$fld=$index->[0][4];
	for $i ( 1 .. $#$index) { 
	  if ($key eq $index->[$i][2]) {
		$fld.=",$index->[$i][4]";
	  } else {
		$sql="create index $key on $ttbl($fld)";
		$fld=$index->[$i][4];
		$key=$index->[$i][2];
		$dbh->do($sql) or print "ERROR:$sql\n";
	  }
	} #for $i

    if ($fld ne '') {
		$sql="create index $key on $ttbl($fld)";
		$fld='';
		$dbh->do($sql) or return;
    }
}


#############################################################################
#  tblExists($dbh,$tbl):0/1
#  useage: tblExists($dbh,'aa')
#  说明:判断table是否存在
#############################################################################
sub tblExists {
  my($dbh,$tbl)=@_;
  my($sql,$sth,$tblnames);
#如果tbl是以xx.yy,则取xx为db
my $db='';
if (index($tbl,'.')!=-1) {
  $db=substr($tbl,0,index($tbl,'.'));
}
$tbl=~s/.*\.//; #去掉db的部分,因为有时传进来的是db.table
$tbl=lc($tbl);
$sql="show tables" if $db eq '';
$sql="show tables in $db" if $db ne '';
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute or die;
$tblnames = $sth->fetchall_arrayref or die ; 

my ($i);
  for $i ( 0 .. $#$tblnames) { 
    return(1) if (lc($tblnames->[$i][0]) eq $tbl);
  } #for $i
return(0);
}

#############################################################################
#  connectDB($conf,$quiet,('chromosome','mpss')):($dbh,$chr,$mpss)
#  useage: ($dbh,$chr)=connectDB('xx.xml',1)
#  就算没有$chr，也需要用 ($dbh)=connectDB(..) ()一定要存在
#  !!不需要事先判断 -e $conf..
#  说明:连接DB 
#  (为统一管理conf.xml文件)
#  !!2011/3/19 若当前文件夹下不存在conf,则去E:\sys\project\_code\XML_conf\下搜
#  #增加输出特定的表名
#  my($dbh,$chrtbl,$mpsstbl)=connectDB($conf,1,('chromosome','mpss')); 或 ($dbh)=connectDB($conf,1);
#############################################################################
sub connectDB{
    my($conf)=shift;
	if (! -e $conf) {
	  $conf=$XMLDIR.$conf;
	  if (! -e $conf) {
		die "$conf not exists!";
	  }
	}
	my($q)=shift;
	my(@tbl)=@_;
	my ($dbh);
	my $xml     = new XML::Simple;
	my $data    = eval {$xml->XMLin($conf , ForceArray => 0)};
	my $db_host= $data->{dbhost};
	my $db_name= $data->{dbname};
	my $db_user_name= $data->{user};
	my $db_password= $data->{password};
	if ($#tbl>=0) {
	  for my $i(0..$#tbl) {
		$tbl[$i]=$data->{$tbl[$i]};
	  }
	}
	$dbh=DBI->connect("DBI:mysql:database=$db_name;host=$db_host",$db_user_name,$db_password);
	if (!$q) {
		if (!$dbh) {
		   die('Database $db_name connect problem!');    
		} else {
		   print "Database $db_name connected!\n";
		}
	}
	return ($dbh,@tbl);
}

#############################################################################
#  @values=getXMLItems($xml,('chromosome','dbname'))
#  说明：得到XML相应项的值
#############################################################################
sub getXMLItems{
    my($xmlfile,@items)=@_;
	if (! -e $xmlfile) {
	  $xmlfile='E:/sys/code/XML_conf/'.$xmlfile;
	  if (! -e $xmlfile) {
		die "xmlfile=$xmlfile not exist!"
	  }
	}
	my $xml     = new XML::Simple;
	my $data    = eval {$xml->XMLin($xmlfile , ForceArray => 0)};
	return () if ($#items<0);
    for my $i(0..$#items) {
	  $items[$i]=$data->{$items[$i]};
    }
	return (@items);
}


#############################################################################
#  getFldsIdx($dbh,$tbl,@flds):@fldsIdx
#  useage: ($a,$b,$notall)=getFldsIdx($dbh,'t_aa',('aa','bb'));
#  说明:返回表中指定字段的下标
#############################################################################
sub getFldsIdx {
	my($dbh,$tbl,@aflds)=@_;
	my($sth,$i,$j,@ids,$notall);
	$sth = $dbh->prepare("desc $tbl") or die $dbh->errstr;
	$sth->execute or die ;
	my $flds= $sth->fetchall_arrayref();
	for $j(0..$#aflds) {
	  $ids[$j]=-1;
      for $i ( 0 .. $#{$flds} ) {	  
		  if (lc($flds->[$i][0]) eq lc($aflds[$j])) {
			  $ids[$j]=$i;
			  last;
		  }
	  }
	}
	#若有一个fld为-1,则返回的数组最后跟一个-1,否则是空
	$notall=0;
	for $i(0..$#ids) {
	  if ($ids[$i]==-1) {
		$notall=-1
	  }
	}
	return $ids[0] if $#ids==0;
	push(@ids,$notall); #添加最后一项
	return @ids;
}

#############################################################################
#  getTblFlds($dbh,$tbl):@fldnames
#  useage: @flds=getTblFlds($dbh,'t_aa');
#  说明:返回表中所有字段名
#############################################################################
sub getTblFlds {
	my($dbh,$tbl)=@_;
	my($sql,$sth,$i,@flds,$tbldefs);
	$sql="desc $tbl";
	$sth = $dbh->prepare($sql) or die $dbh->errstr;
	$sth->execute or die "can't execute the query: $sth->errstr\n";
	$tbldefs = $sth->fetchall_arrayref or die ; 
	for $i ( 0 .. $#{$tbldefs} ) {
		push(@flds,$tbldefs->[$i][0]);
	} 
	return @flds;
}

#############################################################################
#  getFldValues($dbh,$sql,$ncol):@fldValue
#  useage: @fldValue=getFldValues($dbh,'select distinct(chr) from xx order by chr',0);
#  说明:取得某列的值,返回数组. ncol from 0.
#############################################################################
sub getFldValues {
  my($dbh,$sql,$ncol)=@_;
  my($sth,$vs,$i,@values);
  $sth = $dbh->prepare($sql) or die $dbh->errstr;
  $sth->execute or die ;
  $vs = $sth->fetchall_arrayref or die ; 
  for $i(0..$#$vs) {
	push(@values,$vs->[$i][$ncol]);
  }
  return @values;
}

#############################################################################
#  loadFile2Tbl($dbh,$tbl,$file,$ignoreLine):int;
#  useage: $rv=loadFile2Tbl($dbh,$tbl,$file,0)
#  说明:文件导入数据表
#############################################################################
sub loadFile2Tbl {
  my($dbh,$tbl,$file,$ignoreLine)=@_;
  my($sth,$sql,$rv);
  #自动判断 \r\n or \n
  my $mode=fileMode($file);
  my $rn='\r\n'; #win
  if ($mode eq 'unix') {
    $rn='\n';
  }
  my $t='\t';
  if ($ignoreLine<=0) { #\r\n-->\n
    $sql="load data infile '$file' into table $tbl fields terminated by '$t' enclosed by '' Lines Terminated By '$rn'";
  } else {
    $sql="load data infile '$file' into table $tbl fields terminated by '$t' enclosed by '' Lines Terminated By '$rn' ignore $ignoreLine lines";
  }
#  print "$sql\n";
  $sth = $dbh->prepare($sql);
  $rv=$sth->execute or die ;
  return $rv;
}

#############################################################################
#  dotStr2sqlStr($dotstr,$sep=','):sqlstr;
#  useage: $sqlstr=dotStr2sqlStr('chr1,chr2,chr3')
#  说明: 将字符串转换为sql语句的str，输出 'chr1','chr2','chr3'
#############################################################################
sub dotStr2sqlStr {
  my($dotstr,$sep)=@_;
  $sep=',' if !defined($sep) or $sep eq '';
  $dotstr=~s/\s+//g;
  my @dots=split(/$sep/,$dotstr);
  if ($#dots==0) {
	return("'".$dotstr."'");
  } 
  $dotstr=join("','",@dots);
  return("'".$dotstr."'");
}


#****************************************************************************
# Sqlite
#****************************************************************************
#############################################################################
#  getTblFlds_lite($dbh,$tbl):@fldnames
#  useage: @flds=getTblFlds_lite($dbh,'t_aa');
#  说明:返回表中所有字段名
#############################################################################
sub getTblFlds_lite {
	my($dbh,$tbl)=@_;
	my $sth = $dbh->prepare("select sql from sqlite_master where type = \'table\' and tbl_name =\'$tbl\'") or die $dbh->errstr;
	$sth->execute or die ;
	my $flds= $sth->fetchall_arrayref;
	if ($#{$flds}!=0) {
      return ();
	}
    #得到创建表的sql语句,如CREATE TABLE aa (cc int null, bb int null, dd varchar(1024) null)
	my @f;
	my $str=$flds->[0][0];
	my $i1=index($str,'(');
	my $i2=rindex($str,')');
	$str=substr($str,$i1+1,$i2-$i1-1);
	my @ss=split(/,/,$str);
	for my $i(0..$#ss) {
	  $ss[$i]=~s/^\s//;
	  my @ss2=split(/\s/,$ss[$i]);
	  #print "@ss2\n";
	  push(@f,$ss2[0]);
	}
	return @f;
}


#############################################################################
#  getFldsIdx_lite($dbh,$tbl,@flds):@fldsIdx
#  useage: ($a,$b,$notall)=getFldsIdx_lite($dbh,'t_aa',('aa','bb'));
#  说明:返回表中指定字段的下标
#############################################################################
sub getFldsIdx_lite {
	my($dbh,$tbl,@aflds)=@_;
	my($sth,$i,$j,@ids,$notall);
    my @f=getTblFlds_lite($dbh,$tbl);
	if (scalar(@f)<=0) {
      return (-1)x(scalar(@aflds)+1);
	}
	#print "@f\n";
	for $j(0..$#aflds) {
	  $ids[$j]=-1;
      for $i ( 0 .. $#f ) {	  
		  if (lc($f[$i]) eq lc($aflds[$j])) {
			  $ids[$j]=$i;
			  last;
		  }
	  }
	}
	#若有一个fld为-1,则返回的数组最后跟一个-1,否则是空
	$notall=0;
	for $i(0..$#ids) {
	  if ($ids[$i]==-1) {
		$notall=-1
	  }
	}
	return $ids[0] if $#ids==0;
	push(@ids,$notall); #添加最后一项
	return @ids;
}

#############################################################################
#  file2LiteTbl($dbh,$tbl,$file,$ignoreLine):int;
#  useage: $rv=file2LiteTbl($dbh,$tbl,$file,0)
#  说明:文件导入Lite数据表
#############################################################################
sub file2LiteTbl {
  my ($dbh,$tbl,$file,$ignoreLine)=@_;
  my $rv;
  return(0) if !defined($file);
  open (ININ, "<$file" ) or die "Can't open $file!";	
  while ($ignoreLine>0) {
	<ININ>;
    $ignoreLine--;
  }
  my (@items,$l,$sql);
  while (<ININ>) {
    chomp;
	chomp;
	next if $_ eq '';
	@items = split /\t/;
	$l=join (qq{','}, @items);
	$l="\'".$l."\'";
	$sql = "insert into $tbl values($l)";
	$dbh->do( $sql );
	$rv++;
  }  
  $dbh->commit();
  close(ININ);
  return $rv;
}

#############################################################################
#  sql2file($dbh,$sql,$file,$append):rv
#  useage: $rv=sql2file($dbh,$sql,$file,1,'\N')
#  $blankstr用于替代空字段,比如替换成\N,则mysql可以正确导入.
#############################################################################
sub sql2file {
 my($dbh,$sql,$file,$append,$blankstr)=@_;
 my $stm = $dbh->prepare($sql);
 $stm->execute(); 
  if (!$append) {
	open(OO, ">$file") || return 0;  
  } else{
    open(OO, ">>$file") || return 0;  
  } 
my $rv=0;
while (my @row_ary = $stm->fetchrow_array)
{ 
  for my $i(0..($#row_ary-1)) {
	if (!defined($row_ary[$i])) {
	  print OO $blankstr."\t";
	} else {
	  print OO "$row_ary[$i]\t";
	}
  }
  if (!defined($row_ary[$#row_ary])) {
	print OO $blankstr."\n";
  } else {
    print OO "$row_ary[$#row_ary]\n";
  }
  $rv++;
}
return $rv;
}

#############################################################################
##getSqlFlds_Lite($dbh,$tbl,$tblalia,$fldsuf,@excludefields)
##取得类似 a.chr chr_1,a.strand strand_1的字符串
##exclude包含不想要的字段
##Ex. getSqlFlds_Lite($dbh,$tbl,'a','_1',@exclude)
#############################################################################
sub getSqlFlds_Lite {
  my($dbh,$tbl,$prefix,$suf,@ex)=@_;
  $prefix=$tbl if !$prefix;
  $suf='' if !$suf;
  my @flds=getTblFlds_lite($dbh,$tbl);
  my $str='';
  foreach my $fld (@flds) {
	my $no=0;
	if ($#ex!=-1) {
		foreach my $nofld (@ex) {
		  if (lc($fld) eq lc($nofld)) {
			$no=1;
			last;
		  }
		}
	}
    $str.="$prefix.$fld $fld$suf," if !$no;
  }
  $str=substr($str,0,length($str)-1) if $str ne '';
  return($str);
}

#############################################################################
#  tblExists_Lite($dbh,$tbl):0/1
#  useage: tblExists($dbh,'aa')
#  说明:判断table是否存在
#############################################################################
sub tblExists_Lite {
  my($dbh,$tbl)=@_;
  my($sql,$sth,$tblnames);
$tbl=~s/.*\.//; #去掉db的部分,因为有时传进来的是db.table
$tbl=lc($tbl);
$sql="SELECT name FROM sqlite_master WHERE type = \"table\"";
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute or die;
$tblnames = $sth->fetchall_arrayref or die ; 

my ($i);
  for $i ( 0 .. $#$tblnames) { 
    return(1) if (lc($tblnames->[$i][0]) eq $tbl);
  } #for $i
return(0);
}

#############################################################################
#  getTblNames_Lite($dbh):@tables
#  useage: getTblNames_Lite($dbh)
#  说明:取得所有表
#############################################################################
sub getTblNames_Lite {
my($dbh)=shift;
my($sql,$sth,$tblnames);
$sql="SELECT name FROM sqlite_master WHERE type = \"table\"";
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute or die;
$tblnames = $sth->fetchall_arrayref or die ; 

my ($i);
my @ret=();
for $i ( 0 .. $#$tblnames) { 
  push(@ret,$tblnames->[$i][0]) 
} #for $i
return(@ret);
}


#****************************************************************************
# PAT-PA-PAC
#****************************************************************************

#############################################################################
##createPAtbl($dbh,$tbl,$withGFF=0(default)/1,smps=A:B或不提供)
##建PA表
##Ex. createPAtbl($dbh,$tbl,0);
#############################################################################
sub createPAtbl {
  my($dbh,$tbl,$withGFF,$smp)=@_;
  my @smps='tagnum';
  if (defined($smp)) {
	@smps=split(/:/,$smp);	 
  }
  my $txt='';
  for my $s(@smps) {
	$txt=$txt."$s int,";
  }
  $txt=substr($txt,0,length($txt)-1);
  if ($withGFF) {
    $dbh->do("drop table if exists $tbl") or die;
    $dbh->do("create table $tbl(chr varchar(20),strand char(1),coord int,$txt,gff_id int)") or die; 
  } else {
    $dbh->do("drop table if exists $tbl") or die;
    $dbh->do("create table $tbl(chr varchar(20),strand char(1),coord int,$txt)") or die; 
  }
}


#############################################################################
##geneFromGff($dbh,$gfftbl,$otbl,$codon=1/0)
##从gff中提取出gene的start和end
##如果codon=0,则取最开始和最末端
##如果codon=1,则取max(CDS)和min(CDS)
##Ex. $rv=geneFromGff($dbh,'t_gff9','t_gff9_codon_genes',1);
#############################################################################
sub geneFromGff {
  my($dbh,$gfftbl,$otbl,$codon)=@_;
  my $sql;
  if ($codon==0) {
	$sql="select chr,strand,ftr,gene,ftr_start,ftr_end from $gfftbl where ftr not like \'intergenic%\' order by gene";
  } else {
	$sql="select chr,strand,ftr,gene,ftr_start,ftr_end from $gfftbl where ftr=\'CDS\' order by gene";
  }  
  #注意: codon=1的基因数会少于codon=0的,因为codon=1是全部基因,包括psudo的这种,比如27373:33518
  my ($gff,$rv)=execSql($dbh,$sql);
  my $idxs=getIntervals($gff,3);
  my ($startIdx,$endIdx);
  my $tmpfile=getTmpPath(1)."geneFromGff.tmp";
  my $ret=[];
  for my $i(0..$#$idxs) {
    $startIdx=$idxs->[$i][0];
    $endIdx=$idxs->[$i][1];
	my $min=$gff->[$startIdx][4];
	my $max=$gff->[$startIdx][5];
    for my $j(($startIdx+1)..$endIdx) {
	  $min=$gff->[$j][4] if $min>$gff->[$j][4];
	  $max=$gff->[$j][5] if $max<$gff->[$j][5];
	}
   push(@{$ret},[@{$gff->[$startIdx]}]);
   $ret->[$#$ret][4]=$min;
   $ret->[$#$ret][5]=$max;
   $ret->[$#$ret][2]=$ret->[$#$ret][3]
  }	
  saveMtx2File($ret,$tmpfile);
  $ret=[];
  $dbh->do("drop table if exists $otbl") or die;
  $dbh->do("create table $otbl select chr,strand,ftr,gene,ftr_start,ftr_end from $gfftbl where 1<>1") or die;
  $rv=loadFile2Tbl($dbh,$otbl,$tmpfile,0);
  #print "$rv genes by geneFromGff()\n";
  unlink $tmpfile; 
  return($rv);
}


#############################################################################
##cdsLenFromGff($dbh,$gfftbl,$otbl)
##从gff中提取出gene的CDS，计算每个gene的所有CDS的累加长度
#otbl=<gene len>
##Ex. $rv=cdsLenFromGff($dbh,'t_gff9','t_gff9_glen');
#############################################################################
sub cdsLenFromGff {
  my($dbh,$gfftbl,$otbl)=@_;
  my $sql;
  $sql="select chr,strand,ftr,gene,ftr_start,ftr_end from $gfftbl where ftr=\'CDS\' order by gene";
  my ($gff,$rv)=execSql($dbh,$sql);
  my $idxs=getIntervals($gff,3);
  my ($startIdx,$endIdx);
  my $tmpfile=getTmpPath(1)."cdsFromGff.tmp";
  my $ret=[];
  for my $i(0..$#$idxs) {
	my $len=0;
    $startIdx=$idxs->[$i][0];
    $endIdx=$idxs->[$i][1];
	for my $j($startIdx..$endIdx) {
	    $len+=($gff->[$j][5]-$gff->[$j][4]+1);
    }
   push(@{$ret},[$gff->[$startIdx][3],$len]);
  }	
  saveMtx2File($ret,$tmpfile);
  $ret=[];
  $dbh->do("drop table if exists $otbl") or die;
  $dbh->do("create table $otbl select gene,ftr_end len from $gfftbl where 1<>1") or die;
  $rv=loadFile2Tbl($dbh,$otbl,$tmpfile,0);
  #print "$rv genes by geneFromGff()\n";
  unlink $tmpfile; 
  return($rv);
}

#############################################################################
#  cigar2MS
#  useage: cigar2MS(cigar=60M2I5D4M11S includeDI=0)       
#  说明: 由cigar串，得到match和mis的个数
#cigar=60M2I5D4M11S includeDI=0
#M=64,S=11
#cigar=60M2I5D4M11S includeDI=1
#M=71,S=11
#############################################################################
sub cigar2MS {
  my $cigar=shift;
  my $includeDI=shift;
  my ($m,$s)=(0,0);
  my @nums=split(/[A-Za-z]/,$cigar);
  $cigar=~s/^[0-9]+//;
  my @marks=split(/[0-9]+/,$cigar);
  #print "@marks\n@nums\n";
  for my $i(0..$#marks) {
	if ($marks[$i] eq 'M') {
	  $m+=$nums[$i];
	} elsif ($marks[$i] eq 'S') {
	  $s+=$nums[$i];
	} elsif ($includeDI and ($marks[$i] eq 'D' or $marks[$i] eq 'I')) {
	  $m+=$nums[$i];
	}
  }
  return($m,$s);
}

#############################################################################
#  2016/5/24
#  gff文件处理函数
#############################################################################
##读取gff文件全部的chr, ftr和biotype列，得到可能的值
sub getFtrAndBio {
  my $f=shift;
  open(GFF,"<$f") or die "cannot open $f";
  my @ftrs=();
  my @bios=();
  my @chrs=();
  while(my $line=<GFF>) {
	  $line=trim($line);
	  next if $line eq '';
	  next if $line=~/^#/;
	  my @items=split(/\t/,$line);
	  my $chr=$items[0];
	  if ($chr=~/^[1-9]/) {
		$chr="Chr$chr";
	  }
	  my $ftr=$items[2];
	  my $long=$items[8];  
	  my %notes=str2hash($long);
	  my $bio=$notes{'BIOTYPE'};
	  if (!$bio) {
		$bio=$notes{'GENE_BIOTYPE'}; #genebiotype "protein_coding" (mm10.gtf)
	  }
	  push(@ftrs,$ftr);
	  push(@bios,$bio);
	  push(@chrs,$chr);	  
  }
  close(GFF);
  @ftrs=uniq(@ftrs);
  @bios=uniq(@bios);
  @chrs=uniq(@chrs);
  open(OFILE,">${f}.ftr_bio_chr.txt");
  print OFILE "****ftr****\n";
  print OFILE join("\n",@ftrs);
  print OFILE "\n****bio_type****\n";
  print OFILE join("\n",@bios);
  print OFILE "\n****chr****\n";
  print OFILE join("\n",@chrs);
  close(OFILE);
}


#require 需要返回>0的值
1;



