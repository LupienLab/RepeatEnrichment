consensusfiles=$1
peaksdir=$2
peaksext=$4
tedir=$3
teext=$5
outdir=$6
outname=$7


# The first binary matrix is assessing which peak belongs to which sample
## The second binary matrix is assessing which peak overlap with which TE subfamily
for f in $consensusfiles ;
do
    name=$(echo $f | awk -F"/" '{ print $(NF-1) }' )

    echo "----------------"
    echo $name"Binarymat"
    echo "----------------"

    ## For peaks
    ## The second path leads to a directory containing all and only the bedfiles of all the samples you want to run Chromvar for (common end in "_UNIQUE-DNase.component.chrom.start.end.bed" in this case)
    Rscript scripts/01b_createbinarymat_V2.R \
    $f \
    $peaksdir \
    $peaksext \
    $outdir \
    $outname.consensus.Binarymat

    echo "----------------"
    echo $name"Repeat Binarymat"
    echo "----------------"

    ## For repeats
    ## The second path leads to a directory containing all and only the bedfiles of the TE subfamilies (common end in ".bed.sorted")
    #/mnt/work1/software/R/3.4.1/Rscript scripts/createbinarymat_V2.R \
    Rscript scripts/01b_createbinarymat_V2.R \
    $f \
    $tedir \
    $teext \
    $outdir \
    $outname.consensus.Binarymat.repeats

done

# The two binary matrices are generated in the directory where your Consensus set of peaks is located
# Once generated I generally move them to a new directory which contains a subdirectory called "chromvar". The chromvar directory is necessary for the following step
