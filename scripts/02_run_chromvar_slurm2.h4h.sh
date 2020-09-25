peaksbinarymat=$1
peaksbinarymatrepeats=$2
outname=$3
outdir=$4

Rscript scripts/02b_chromvar_h4h.R $peaksbinarymat $peaksbinarymatrepeats $outname $outdir
