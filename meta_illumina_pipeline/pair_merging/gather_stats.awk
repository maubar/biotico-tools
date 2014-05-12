BEGIN{ ok = 0; reads = 0; no_algn=0; lowq=0;bad_r=0}
$3 == "OK"{ ok += $4}
$3 == "READS"{ reads += $4 }
$3 == "LOWQ"{ lowq += $4 }
$3 == "NOALGN"{ no_algn += $4 }
$3 == "BADR"{ bad_r += $4 }
END{
	print "READS:\t",reads;
	print "OK:\t",ok;
	print "LOWQ:\t",lowq;
	print "NOALGN:\t",no_algn;
	print "BADR:\t",bad_r;
}
