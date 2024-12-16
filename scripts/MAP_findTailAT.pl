#!/usr/bin/perl -w


use strict;
use Getopt::Long;
use File::Basename;


##########################################################################################
#  parameters & usage
##########################################################################################
my %opts=();
GetOptions(\%opts,"h","in=s","poly=s","ml=i","mp=i","mg=i","mm=i","mtail:i","mper:s","mr:i","odir:s","suf:s","reg:i","deep:s","oraw:s","debug:s","bar:s","review:s");
my $USAGE=<< "_USAGE_";
Usage: 
  0) Test different parameters -debug T open
  -debug T -oraw T

  1) Recommended parameters for poly(T) stretch type of 3'seq dataset
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly T -ml 25 -mp 6 -mg 5 -mm 2 -mr 2 -mper 0.75 -mtail 8 -deep T -reg 1 -odir "f:/" -suf "" 
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly T -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "" -reg 1
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly "A|T" -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "AT" -reg 1

  Do not output the original sequence
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly "A&T" -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "A&T" -reg 1 -oraw F

  debug output various situations
  MAP_findTailAT.pl -in E:/sys/code/testdata/arab1.fastq -poly "A&T" -ml 25 -mp 8 -mg 8 -mm 2 -odir "f:/" -suf "A&T" -reg 1 -oraw F -debug T

-h=help
-in=input fa or fq file
-poly=A/T/A&T/A|T
  A|T If there are both A and T, when tail_lenA>=lenT, keep the result of A; a sequence can only contain A or T, not both
  A&T Find A and T at the same time, the two are irrelevant, a sequence may have both A and T
-ml=min length after trim
-mg=margin from the start (poly=T) or to the end (poly=A) (=5)
-mm=mismatch between TxxTTT (=2)
-mr=minT in reg (=3)
-mp=min length of succesive poly (=8)
-mtail=min length of trimmed tail (=8)
-mper=min percent of A/T in trimmed tail (=0.75)
-deep=T(default)/F If is T tail, then if reg does not match, search deeply, such as this typeTTTTTCTTTTCTCTTTTTTTT
-odir=output path (default is the same as input)
-suf=suffix, default: xx.suf.T/A.fq
-reg=1(loose, default)/2
  1='^.{0,mg}?(T{mr,}[^T]{0,mm})*T{mp,}([^T]{0,mm}T{mr,})*'
  2='^.{0,mg}?(T{mp,}([^T]{0,mm}T{mr,})*)'
-oraw=T(default)/F; output raw reads for trimmed A/T file (.A.raw.fq)
-bar=Last bast position of three barcode 
-review=T (default) or F, if T then when no match was found for reg or deep, do deeper search for tails like AACCCCTTTTTTTTTTTTTTTTT...
_USAGE_


#############################################################################
#  invoke the program                               
#############################################################################
die $USAGE if $opts{'h'}||!$opts{'in'}||!$opts{'poly'}||!$opts{'ml'}||!$opts{'mp'}||!$opts{'mg'}||!$opts{'mm'};

my $poly=$opts{'poly'};
my $in=$opts{'in'};
my $ml=$opts{'ml'};
my $mp=$opts{'mp'};
my $mg=$opts{'mg'};
my $mm=$opts{'mm'};
my $mr=$opts{'mr'};
my $mtail=$opts{'mtail'};
my $mper=$opts{'mper'};
my $deep=$opts{'deep'};
my $debug=$opts{'debug'};
my $suf=$opts{'suf'};
my $odir=$opts{'odir'};
my $reg=$opts{'reg'};
my $oraw=$opts{'oraw'};
my $bar=$opts{'bar'};
my $review=$opts{'review'};
$reg=1 if !$reg;
die "reg=1/2" if $reg!=1 and $reg!=2;


if (!$mtail) {
  $mtail=8;
}
if (!$mper) {
  $mper=0.75;
}
if (!$mr) {
  $mr=3;
}
if (!$deep | $deep eq 'T') {
  $deep=1;
} else {
  $deep=0;
}
if (!$review | $review eq 'T') {
  $review=1;
} else {
  $review=0;
}
if (!$debug | $debug ne 'T') {
  $debug=0;
} else {
  $debug=1;
}
if(!$bar){
 $bar=0;	
}
$oraw=1 if !$oraw or $oraw ne 'F';
$oraw=0 if $oraw eq 'F';

print "poly=$poly\ninfile=$in\noutput dir=$odir\nmin length after trim(ml)=$ml\nmin continue tail length(mp)=$mp\nmax margin to the end(mg)=$mg\nmax mismatches between xTTT(mm)=$mm\n";
print "min reg T(mr)=$mr\nmin tail length(mtail)=$mtail\nmin T/A percent(mper)=$mper\n";
print "deep=$deep\n";
print "DEBUG MODE, output to .test.file too \n" if $debug;

die "error poly(=A T A|T A&T)" if $poly ne 'A' and $poly ne 'T' and $poly ne 'A|T' and $poly ne 'A&T';
my ($findT,$findA)=(0,0);
$findA=1 if ($poly=~ m/A/);
$findT=1 if ($poly=~ m/T/);
#print "$findT\n";

if ($deep & $findA) {
  die "TODO: The deep and findA procedures are not yet implemented...";
}

die "only support fastq file" if seqFormat($in) ne 'fq';
open my $fh1, "<$in" or die "cannot read input file ($in)\n";
#open IN, "<$in" or die "cannot read input file ($in)\n";
$suf=".$suf" if $suf;

my ($cntNotailT,$cntShortT,$cntTotal,$cntFinalT,$cntNotailA,$cntShortA,$cntFinalA,$cntMissA,$cntMissT,$cntBadT,$cntBadA)=(0,0,0,0,0,0,0,0,0,0,0);

my $inpath=getDir($in,0); # x/y/
my $infname=getFileName($in,1); #
if (!$odir) {
  $odir=$inpath;
}

if ($findA) {
  open OA,">$odir$infname$suf.A.fq";
  
 if($debug){
  	open TESTA,">$odir$infname$suf.testA.fq";
  	}
#  if($oraw){
#  	open RAWA,">$odir$infname$suf.A.raw.fq" ;
#  	}
}


if ($findT) {
  open OT,">$odir$infname$suf.T.fq";
  if($debug){
#  	open TESTT,">$odir$infname$suf.T.DEBUG.fq";
#  	open BEDT,">$odir$infname$suf.badT.fq";
#  	open NOT,">$odir$infname$suf.noTail.fq";
 	open SHORTT,">$odir$infname$suf.shortT.fq";
  	}
#  if($oraw){
#  	open RAWT,">$odir$infname$suf.T.raw.fq";
#  	}
}

#print "###############$findT\t$debug\t$oraw\n";

my $regT="^.{0,$mg}?(T{$mr,}[^T]{0,$mm})*T{$mp,}([^T]{0,$mm}T{$mr,})*";
if ($reg==2) {
  $regT="^.{0,$mg}?T{$mp,}([^T]{0,$mm}T{$mr,})*";
}
print "regT=$regT\n";
#my $regA="(A{$mr,}[^A]{0,$mm})*A{$mp,}([^A]{0,$mm}A{$mr,})*.{0,$mg}?\$";
my $regA="(A{$mr,}[^A]{0,$mm})*A{$mp,}[^A]{0,1}A+([^A]{0,$mm}A{$mr,})*.{0,$mg}?\$";

if ($reg==2) {
  $regA="A{$mp,}([^A]{0,$mm}A{$mr,})*.{0,$mg}?\$";
}
print "regA=$regA\n";

if ($deep & $findT) {
  print "deep reg=^.{0,$mg}?(T{1,}[^T]{0,2})*T{8,}\n";
}
if($review & $findT){
   print "review_reg=T{$mp,}([^T]{0,$mm}T{$mr,}){0,}\n";	
}

my ($md1,$md2,$md3)=(0,0,0);


my @fqkeys = qw/id seq pl qual/;
while(!eof $fh1) {
	my (%e1);
	@e1{@fqkeys} = readfq($fh1);
	my $tempID = $e1{'id'};
	my @spl = split(/\s+/,$tempID);
	$e1{'id'} = $spl[0];

	$cntTotal++;
	my ($seqT,$seqA,$qualT,$qualA,$trimS)=('','','','','');
	my ($haveT,$haveA,$shortT,$notailT,$shortA,$notailA,$badTailT,$badTailA)=(0,0,0,0,0,0,0,0);
	if ($findT) {
	    $seqT=$e1{'seq'};
 		if ($seqT=~s/$regT//) {
		  my $length=length($seqT);
		  #my $Ts=substr($e1{'seq'},0,length($e1{'seq'})-$length); 
		  my $Ts=substr($e1{'seq'},$bar,length($e1{'seq'})-$length-$bar);
		  #print "$Ts\n";#Seven test
		  #my $Tcnt = $Ts =~ tr/T/T/; 
		  my $Tcnt = $Ts=~ tr/T/T/;
		  ##print "Yes: $Ts.$Tcnt\n";
		  if($length<$ml){
				$shortT=1;
				$badTailT=0;
				#print "model1\t$seqT\t$Ts\n";
				 #print SHORTT "$e1{'seq'}\n" if $debug; #seven
			#print TESTT $Ts."***SHR***".$seqT."\n"  if $debug;	  	
		  } elsif ((length($e1{'seq'})-$length-$bar)<$mtail or $Tcnt/length($Ts)<$mper) { #tail过短 or T%<mper
			$badTailT=1;
		  $shortT=0;
			#print TESTT $Ts."***BAD***".$seqT."\n"  if $debug;
		  } else {
		  	++$md1;
		  	$haveT=1;
		  	$trimS=$Ts; #trimmed part sequence
		    $qualT=substr($e1{'qual'},length($e1{'qual'})-$length,$length);
		     print OT "$e1{'id'}\_$trimS\n$seqT\n$e1{'pl'}\n$qualT\n";
		     $cntFinalT++;
				 $badTailT=0;
				 $shortT=0;			     
		  }
		} 

    if (!$haveT & $deep) { 
          $seqT=$e1{'seq'};
		      my $deepreg="^.{0,$mg}?(T{1,}[^T]{0,2})*T{$mp,}";
		      if ($seqT=~s/$deepreg//){
			   		 my $length=length($seqT);
			 			 #my $Ts=substr($e1{'seq'},0,length($e1{'seq'})-$length); 
			 			 my $Ts=substr($e1{'seq'},$bar,length($e1{'seq'})-$length-$bar);
			 			   #print "$Ts\n";
			  	  my $Tcnt = $Ts =~ tr/T/T/; 
		      		##print "Deep: $Ts.$Tcnt\n";
		      		
		      		
		     if   ($length<$ml) { 
			  			$shortT=1;
			  			$badTailT=0;
			  				#print "model2\t$seqT\t$Ts\n";
			  			# print SHORTT "$e1{'seq'}\n" if $debug; #seven	`
			    } elsif	((length($e1{'seq'})-$length-$bar)<$mtail or $Tcnt/length($Ts)<$mper) {
			  			$badTailT=1;
			  			$shortT=0;
			  	}else{
			  		  ++$md2;
				      $haveT=1;
				      $trimS=$Ts; #trimmed part sequence
				      $qualT=substr($e1{'qual'},length($e1{'qual'})-$length,$length);
				      print OT "$e1{'id'}\_$trimS\n$seqT\n$e1{'pl'}\n$qualT\n";
				      $cntFinalT++;
				      #print TESTT $Ts."***DEP***".$seqT."\n"  if $debug;
				      $badTailT=0;
				      $shortT=0;		  	
			  		}	      	 
		     }
	  }
	     
	 if(!$haveT & $review){
	     	   $seqT=$e1{'seq'};
	     	   my $deepmod="^.{0,$mg}?T{$mp,}([^T]{0,$mm}T{$mr,}){0,}";
	     	   if($seqT=~ m/$deepmod/){
	     	   	
	     	     #print "1$seqT\n";
	     	     my $pos_start=$-[0];
	     	     my $match_length =length$&;
             my $pos_end = $pos_start+$match_length;
  
	     	   	 $seqT =substr($seqT,$pos_end,length($seqT)-$pos_end);
	     	   	 #print "2$seqT\n";
	     	   	 my $length=length($seqT);
	     	   	 
	     	   	 my $Ts=substr($e1{'seq'},$bar,$pos_end-$bar);#polyT
	     	   	 #print "3$Ts\n";
	     	   	 my $Tcnt = $Ts =~ tr/T/T/; #
	     	   	 #print "3$Tcnt\n";
	     	  
	     	  if ($length<$ml) { #seq is too short
			  			$shortT=1;
			  			$badTailT=0;
			  					#print "model3\t$seqT\t$Ts\n";
			  			 #print SHORTT "$e1{'seq'}\n" if $debug; #seven
			    } elsif ( length($Ts)<$mtail or $Tcnt/length($Ts)<$mper) {#tail过短 or T%<mper	
			  			$badTailT=1;
			  		  $shortT=0;
			  	} else{
			   	  ++$md3;
	     	   	$haveT=1;
	     	   	$trimS=$Ts; #trimmed part sequence
	     	   	$qualT=substr($e1{'qual'},$pos_end,length($e1{'qual'})-$pos_end);
	     	   	print OT "$e1{'id'}\_$trimS\n$seqT\n$e1{'pl'}\n$qualT\n";
	     	   	$cntFinalT++;
				    $badTailT=0;
				    $shortT=0;
			   	}
	     	   	}
	     }
	     

   if (!$haveT & !$shortT) {
		  $notailT=1;	
		  $cntNotailT++;	
		  print OT "$e1{'id'}\_$trimS\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";  
	 } 
	 if($badTailT) {
	 	   $cntBadT++;
	 }
	 if($shortT) {
			$cntShortT++;
			 print SHORTT "$e1{'id'}\n$e1{'seq'}\n" if $debug; #seven
	  }
	
	 
	 
	}#end find poly(T)


	if ($findA) {
	    $seqA=$e1{'seq'};
		if ($seqA=~s/$regA//) {
	   my $length=length($seqA);
		  #print TESTA "$seqA*******".substr($e1{'seq'},$length,length($e1{'seq'})-$length)."\n" if $debug;
		  my $As=substr($e1{'seq'},$length,length($e1{'seq'})-$length); #截取下来的A片段
		  my $Acnt = $As =~ tr/A/A/; #计算A数
		  
		  if ($length<$ml) { #default ml=20 or 40
		   $shortA=1;
		   $badTailA=0;
			 #$cntShortA++;
			 ##print OA "$e1{'id'}\\_$trimS\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n"; 
		  } elsif (length($e1{'seq'})-$length<$mtail or $Acnt/length($As)<$mper) { #tail过短
			$badTailA=1;
			#$cntBadA++;
		  #print TESTA "Model1-Short:*******$e1{'seq'}\n" if $debug;
		  } else {
			$haveA=1;
			$shortA=0;
			$trimS=$As; #trimmed part sequence
			$qualA=substr($e1{'qual'},0,$length);
			print OA "$e1{'id'}\_$trimS\n$seqA\n$e1{'pl'}\n$qualA\n";
			$cntFinalA++;
		  }
		} else {
		  #print OA "$e1{'id'}\_$trimS\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";  
		  $notailA=1;
		  $badTailA=0;
		  $shortA=0;
		  #$cntNotailA++;
		  #print TESTA "Model1-No:*******$e1{'seq'}\n" if $debug;
		}
		
		if(!$haveA & $review){
			
				 my $deepmod="A{10,}([^A]{0,2}A{3,})?.{0,50}?\$";
			   $seqA=$e1{'seq'};
			 	    if ($seqA=~s/$deepmod//) {
	         	my $length=length($seqA);
		 		 #print TESTA "$seqA*******".substr($e1{'seq'},$length,length($e1{'seq'})-$length)."\n" if $debug;
		  	   	my $As=substr($e1{'seq'},$length,length($e1{'seq'})-$length); #截取下来的A片段
		  	   	my $Acnt = $As =~ tr/A/A/; #计算A数
		  
		  		if ($length<$ml) { #default ml=20 or 40
		   				$shortA=1;
		   				$badTailA=0;
			 				#$cntShortA++;
			 	##print OA "$e1{'id'}\\_$trimS\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n"; 
		  		} elsif (length($e1{'seq'})-$length<$mtail) { #tail过短
						$badTailA=1;
						$shortA=0;
						#$cntBadA++;
		        #print TESTA "Model2-Short:*******$e1{'seq'}\n" if $debug;
		  		} else {
						$haveA=1;
						$trimS=$As; #trimmed part sequence
						$qualA=substr($e1{'qual'},0,$length);
						print OA "$e1{'id'}\_$trimS\n$seqA\n$e1{'pl'}\n$qualA\n";
						$cntFinalA++;
		  		}
		  }else{
		  	 $notailA=1;
			   #print TESTA "Model1-No:*******$e1{'seq'}\n" if $debug;
			}
			
			
		}
	
	
	  if (!$haveA & !$shortA) {
		  $notailA=1;	
		  $cntNotailA++;	
		  print OA "$e1{'id'}\_$trimS\n$e1{'seq'}\n$e1{'pl'}\n$e1{'qual'}\n";  
	   } 
		
		if(!$haveA & $shortA){
			  $cntShortA++;
			  print TESTA "Short:*******$e1{'seq'}\n";
			}
		if(!$haveA & $badTailA & !$shortA){
			$cntBadA++;
			print TESTA "BadA:*******$e1{'seq'}\n";
			}
		if(!$haveA & $notailA & !$shortA){
			$cntNotailA++;
			print TESTA "noTail:*******$e1{'seq'}\n";
		}
		
		
		
		
		
		
		
	} #~findA






}#end while 
close($fh1);

if($findA){
   close(OA);
#   if($debug){
#   	close(TESTA);
#   	}	
#   	if($oraw){
#   	close(RAWA);	
#   	}
}


if($findT){
   close(OT);
  if($debug){
#   	close(TESTT);
#   	close(NOT);
#   	close(BEDT);
   	close(SHORTT);
   	}	
#   	if($oraw){
#   	close(RAWT);	
#   	}
}



print "\ntotal\t$cntTotal\n";
if ($findA) {
  print "[polyA]\nfinal\t$cntFinalA\nnotail\t$cntNotailA\ntooshort\t$cntShortA\nbadTail\t$cntBadA\nmissby(A|T)\t$cntMissA\n";
  if ($cntTotal!=$cntFinalA+$cntNotailA+$cntShortA+$cntMissA+$cntBadA) {
	print "cntTotal!=cntFinalA+cntNotailA+cntShortA+cntMissA+cntBadA\n";
  }
}
if ($findT) {
  print "[polyT]\nfinal\t$cntFinalT\nnotail\t$cntNotailT\ntooshort\t$cntShortT\nbadTail\t$cntBadT\nmissby(A|T)\t$cntMissT\n";
  print "region\t$md1\ndeep\t$md2\nreview\t$md3\n";
  if ($cntTotal!=$cntFinalT+$cntNotailT+$cntShortT+$cntMissT+$cntBadT) {
	print "cntTotal!=cntFinalT+cntNotailT+cntShortT+cntMissT+cntBadT\n";
  }
}


###############################################
sub readfq {
        my $fh = pop @_;
        my @entry = ();
        for (qw/id seq pl qual/) {
                my $line = readline($fh); chomp $line; 
                push @entry, $line; warn "line empty" if !$line;
        }
        return @entry;
}


#############################################################################
#  getDir(afilename,nobar=0/1) 
#  useage: getDir('xx/xx.txt'); getFileName('xx/xx.txt',1)        
#  ËµÃ÷: ·µ»ØÂ·¾¶
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


#############################################################################
#  getFileNames($dir,$pattern,$NOTMATCH):@files
#  useage: @files=getFileNames('.','\.pl'); @files=getFileNames('.','\.pl',1);      
#  ËµÃ÷: ·µ»ØÎÄ¼þÃû°üº¬Â·¾¶,ÊäÈëÂ·¾¶¿Éº¬/»ò²»º¬,Æ¥ÅäÀ©Õ¹ÃûÒªÓÃÈç"\.txt"»ò txt$; 
#       Èô$NOTMATCH·Ç£°,Ôò±íÊ¾²»Æ¥ÅäÄ£Ê½²ÅÊä³öÎÄ¼þ
#  ÈôÒªÆ¥Åä xxtrain1.arff; train222.arff; ÔòÓÃ pat="train.*\.arff"
#  ÈôÒªÆ¥ÅäÒÔtrain¿ªÍ·ÒÔ.arff½áÎ²µÄ;ÔòÓÃ pat="^train.*\.arff$"
#  ÒÔarff½áÎ²µÄ arff$
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


#****************************************************************************
# Sequencing
#****************************************************************************

#############################################################################
#  seqFormat($seqfile):fa/fq/''
#  useage: $format=seqFormat($file)
#  ËµÃ÷: ÅÐ¶ÏÐòÁÐÎÄ¼þÊÇfa»òfq¸ñÊ½,Ö»Í¨¹ýµÚ1ÐÐµÄ±êÌâÊÇ>»¹ÊÇ@À´ÅÐ¶Ï¸ñÊ½
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
#  getFileName(afilename,noExt=0/1) 
#  useage: getFileName('xx/xx.txt'); getFileName('xx/xx.txt',1)        
#  ËµÃ÷: ·µ»ØÎÄ¼þÃû£¨²»º¬Â·¾¶£©
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

