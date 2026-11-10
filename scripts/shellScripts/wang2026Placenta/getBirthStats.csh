#!/bin/tcsh
#SBATCH --job-name=wang26BirthStats
#SBATCH --output=logs/wang26BirthStats.out
#SBATCH --error=logs/wang26BirthStats.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300G

# Leave unmatched globs as-is instead of aborting with "No match."
set nonomatch

# ---------------------------------------------------------------------------
# Paths and settings
# ---------------------------------------------------------------------------
set bedMerge        = $HOME/go/bin/bedMerge
set kent            = /hpc/group/vertgenlab/softwareShared/kent/kent.Jul.30.2021

set workDir         = /work/hcn4/260630_vertCons_wd/thirdTriPlacenta_snATAC/Wang_2026
set rawPeaks        = $workDir/narrowPeaks/raw
set nonOverlap      = $workDir/narrowPeaks/nonOverlapping
set uniqueFile      = $workDir/wang2026Placenta.4colUnique.bed
set workingFile     = $workDir/wang2026Placenta.4colMerged.bed
set outFile         = $workDir/wang2026PlacentaStatsMerged.txt
set birthStatsDir   = $workDir/birthStats

set db              = hg38      # reference
set allSpecies      = /hpc/group/vertgenlab/christi/vertCons/60way/allSpeciesBeds/$db.allSpecies.bed
set dataFile        = $uniqueFile
set speciesNodeFile = /hpc/group/vertgenlab/christi/vertCons/60way/speciesNodeFiles/$db.speciesNodes.list
set threshold       = 0.33
set holepunch       = /hpc/group/vertgenlab/hailey/vertCons/code/vertCons_revision/scripts/shellScripts/wang2026Placenta 

# ---------------------------------------------------------------------------
# 1.0 Remove any peaks that overlap within 200bp
# ---------------------------------------------------------------------------
echo PROCESSING PEAKS...
mkdir -p $workDir $nonOverlap

foreach f ($rawPeaks/*)
    set base = $f:t:r                      # basename minus last extension
    echo $f

    # Collapse each peak to its summit (narrowPeak col 10 = offset from start),
    # then extend 200bp each way (400bp total). Falls back to the midpoint when
    # no summit is given (col 10 = -1 or missing). Skips track/comment lines;
    # start is clamped at 0 near chromosome ends.
    awk 'BEGIN{FS=OFS="\t"} /^(track|#)/{next} {if (NF>=10 && $10>=0) c=$2+$10; else c=int(($2+$3)/2); s=c-200; if(s<0)s=0; $2=s; $3=c+200; print}' \
        $f > $nonOverlap/$base.resized.bed

    $bedMerge -pad 200 $nonOverlap/$base.resized.bed $nonOverlap/$base.bed
    rm -f $nonOverlap/$base.resized.bed
    gzip -f $nonOverlap/$base.bed          # -f: overwrite .gz from a previous run
end

set peakFiles = ($nonOverlap/*.gz)
if (! -e "$peakFiles[1]") then
    echo "ERROR: no .gz files found in $nonOverlap"
    exit 1
endif

echo DONE PROCESSING PEAKS


# ---------------------------------------------------------------------------
# 2.0 Get catted atlas bed file   (was step 3; must run before the stats step)
## From /hpc/group/vertgenlab/christi/vertCons/60way/scriptPipeline/9.nameAndCatZippedAtlas.csh
# ---------------------------------------------------------------------------
echo CATTING ATLAS BED FILES...
rm -f $uniqueFile

foreach f ($peakFiles)
    set name = $f:t:r:r
    # Always write exactly 4 columns: chrom, start, end, name.
    # If the file has no name column (BED3), make one from the row number.
    zcat $f | tr " " "\t" | awk -v n="$name" 'BEGIN{FS=OFS="\t"} {id = (NF>=4) ? $4 : "peak"NR; print $1, $2, $3, id".line"NR"."n}' >> $uniqueFile
end

$bedMerge -keepAllNames $uniqueFile $workingFile

echo DONE CATTING ATLAS BED FILES


# ---------------------------------------------------------------------------
# 3.0 Get atlas stats   (was step 2)
## From /hpc/group/vertgenlab/christi/vertCons/60way/scriptPipeline/8.atlasStats.csh
# ---------------------------------------------------------------------------
echo GETTING ATLAS STATS...
rm -f $outFile

foreach f ($peakFiles)
    set name  = $f:t:r:r
    set lines = `grep -c -- $name $workingFile`

    echo $name
    echo $lines

    printf "%s\t%s\n" $name $lines >> $outFile
end

echo DONE GETTING ATLAS STATS


# ---------------------------------------------------------------------------
# 4.0 Calculate ancestor
## From /hpc/group/vertgenlab/christi/vertCons/60way/scriptPipeline/10.calculateAncestor.csh
# ---------------------------------------------------------------------------
echo CALCULATING ANCESTOR...
cd $workDir
rm -rf coverageCalls.$db
mkdir coverageCalls.$db

# Needs kentutils (overlapSelect). helper.punchHolesInBed.pl must be in $workDir.
$kent/overlapSelect -nonOverlapping $allSpecies $dataFile $db.specific.bed
$kent/overlapSelect -mergeOutput $allSpecies $dataFile stdout \
    | $holepunch/helper.punchHolesInBed.pl /dev/stdin coverageCalls.$db/tmp.a.bed

foreach x (`cat $speciesNodeFile`)
    set species = $x:r
    set node    = $x:e
    echo $species $node

    $kent/overlapSelect -aggregate -overlapThreshold=$threshold \
        /work/hcn4/260630_vertCons_wd/$db.wholeGenomeAlignments/$db.$species.bed \
        coverageCalls.$db/tmp.a.bed coverageCalls.$db/tmp.out.bed

    cat coverageCalls.$db/tmp.out.bed >> coverageCalls.$db/node"$node".bed
    cut -f 4 coverageCalls.$db/tmp.out.bed | awk -v "NODE=$node" '{ print $1"\t"NODE }' >> coverageCalls.$db/cneNode.tsv

    $kent/overlapSelect -nonOverlapping coverageCalls.$db/node"$node".bed \
        coverageCalls.$db/tmp.a.bed coverageCalls.$db/tmp.b.bed

    mv coverageCalls.$db/tmp.b.bed coverageCalls.$db/tmp.a.bed
    rm coverageCalls.$db/tmp.out.bed
end

mv coverageCalls.$db/tmp.a.bed coverageCalls.$db/refSpecific.bed
wc -l coverageCalls.$db/refSpecific.bed

echo DONE CALCULATING ANCESTOR


# ---------------------------------------------------------------------------
# 5.0 Calculate birth node stats
## From /hpc/group/vertgenlab/christi/vertCons/60way/scriptPipeline/11.birthNodeStats.csh
# ---------------------------------------------------------------------------
echo CALCULATING BIRTH NODE STATS...
mkdir -p $birthStatsDir

foreach n (coverageCalls.$db/node*)
    set node       = $n:t:r
    set birthStats = $birthStatsDir/$node.stats.txt

    printf "Node\tCell Type\tRegions\n" > $birthStats

    foreach f ($peakFiles)
        set name  = $f:t:r:r
        set lines = `grep -c -- $name $n`

        echo $node
        echo $name
        echo $lines

        printf "%s\t%s\t%s\n" $node $name $lines >> $birthStats
    end
end

cd $birthStatsDir
cat node*.txt >wang2026Placenta_gains.txt
gawk -i inplace '!($0 == "Node\tCell Type\tRegions" && seen++)' wang2026Placenta_gains.txt 

echo DONE CALCULATING BIRTH STATS


# ---------------------------------------------------------------------------
# 6.0 Calculate death node
## From /hpc/group/vertgenlab/christi/vertCons/60way/scriptPipeline/12.calculateElementDeathAncestors.csh
# ---------------------------------------------------------------------------
echo CALCULATING DEATH NODE...

set deathNodeDir      = $workDir/deathNode/$db
set intervalOverlap   = $HOME/go/bin/intervalOverlap
set alignmentDir      = /work/hcn4/260630_vertCons_wd/$db.wholeGenomeAlignments
set coverage          = $workDir/coverageCalls.$db
set deathSpeciesNodes = /hpc/group/vertgenlab/christi/vertCons/60way/speciesNodeFiles/$db.speciesNode.list
set maxNode           = 17

mkdir -p $deathNodeDir
rm -f $deathNodeDir/*.bed              # clear output from a previous run

# Add a column with the birth node to each node file
foreach f ($coverage/node*.bed)
    set node = $f:t:r
    if ("$node" =~ *.labeled) continue  # skip already-labeled files
    echo $node
    awk -v node="$node" 'BEGIN{OFS="\t"} {print $0, node}' $f > $f:r.labeled.bed
end

# For each birth node (oldest to youngest), find elements with no alignment
# in each species that branched off after that node
foreach i (`seq $maxNode -1 0`)
    set nodeFile = $coverage/node$i.labeled.bed
    if (! -e $nodeFile) continue        # no elements born at this node

    foreach sN (`cat $deathSpeciesNodes`)
        set species = $sN:r
        set node    = $sN:e
        if ("$species" == "$db") continue
        if ($node >= $i) continue

        echo "node: $node  i: $i  species: $species"
        $intervalOverlap -nonOverlap $alignmentDir/$db.$species.bed $nodeFile \
            $deathNodeDir/node$i.$sN.dead.nonOverlap.bed
    end
end

# Add columns with the death node and species
foreach f ($deathNodeDir/node*.*.*.dead.nonOverlap.bed)
    set death   = $f:t:r:r:r:e
    set species = $f:t:r:r:r:r:e
    set prefix  = $f:r:r
    echo $death $species

    awk -v death="node$death" -v species="$species" 'BEGIN{OFS="\t"} {print $0, death, species}' $f > $prefix.allLabels.bed
end

echo DONE CALCULATING DEATH NODE


# ---------------------------------------------------------------------------
# 7.0 Calculate death stats
## From /hpc/group/vertgenlab/christi/vertCons/60way/scriptPipeline/13.deathStats.csh
# ---------------------------------------------------------------------------
echo CALCULATING DEATH STATS...

set deathOut = $workDir/allDeathStats.updated.tsv

# Combine each species' dead elements across birth nodes
foreach sN (`cat $deathSpeciesNodes`)
    set s = $sN:r
    if ("$s" == "$db") continue

    set deadFiles = ($deathNodeDir/node*.$s.*.dead.allLabels.bed)
    if (! -e "$deadFiles[1]") then
        echo "No death files for $s"
        continue
    endif

    cat $deadFiles > $deathNodeDir/$s.dead.bed
    # bedMerge needs plain 4-column BED (the label columns would break it)
    cut -f1-4 $deathNodeDir/$s.dead.bed > $deathNodeDir/$s.dead.4col.bed
    $bedMerge -keepAllNames $deathNodeDir/$s.dead.4col.bed $deathNodeDir/$s.dead.merged.bed
end

# For each cell type and species: elements born at or before the species'
# branch point, and how many of those are dead in that species
rm -f $deathOut

foreach f ($peakFiles)
    set name = $f:t:r:r

    foreach sN (`cat $deathSpeciesNodes`)
        set species = $sN:r
        set node    = $sN:e
        if ("$species" == "petMar3" || "$species" == "$db") continue

        set sum = 0
        foreach i (`seq $maxNode -1 $node`)
            if (-e $coverage/node$i.bed) then
                set tmp = `grep -c -- $name $coverage/node$i.bed`
                @ sum += $tmp
            endif
        end

        set dead = 0
        if (-e $deathNodeDir/$species.dead.bed) then
            set dead = `grep -c -- $name $deathNodeDir/$species.dead.bed`
        endif

        printf "%s\t%s\t%s\t%s\n" $name $species $sum $dead >> $deathOut
    end
end

echo DONE!
