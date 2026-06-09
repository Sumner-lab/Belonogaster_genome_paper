#!/usr/bin/perl

#==============================================================================
# ChopGO.pl
#------------------------------------------------------------------------------
# Gene Ontology (GO) enrichment for one or more gene lists, using topGO in R.
#
# Given a list of genes (optionally split into subgroups via a second column)
# and a GeneID -> GO mapping (either from VASTDB by species code, or a custom
# --GO_file), the script:
#   1. Builds the gene-to-GO background for R (all genes, or a -bg subset).
#   2. Writes per-group topGO scripts (BP, MF and CC ontologies).
#   3. Runs R to perform classic Fisher enrichment with multiple-testing
#      corrections, writing a "<group>_TopGo_results_ALL.tab" table per group.
#   4. Optionally generates histogram plots via MakeHist_fromChartGO-10.pl.
#
# Pre-requisites: R (with the topGO Bioconductor package).
#
# Questions:
#   Chris Wyatt   (cw13722@gmail.com)
#==============================================================================

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Cwd qw(abs_path cwd);
use File::Path qw(make_path);

#------------------------------------------------------------------------------
# Initialise paths
#------------------------------------------------------------------------------
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//;

#------------------------------------------------------------------------------
# Default option values
#------------------------------------------------------------------------------
my $GO_file;
my $input;
my $background;
my $species;
my $helpFlag      = 0;
my $path_to_DB    = "$binPath/DB";
my $plot          = 1;
my $pval_cutoff   = 0.05;
my $num_cutoff    = 10;
my $plot_only     = 0;
my $plot_open     = 0;
my $install_topGO;
my $meth          = "none";
my $sort_enrich   = 0;
my $min_genes     = 2;
my $max_genes     = 100000;
my $infer         = 0;

#------------------------------------------------------------------------------
# Parse command-line options
#------------------------------------------------------------------------------
GetOptions(
    "help"            => \$helpFlag,
    'i=s'             => \$input,
    "bg=s"            => \$background,
    'sp=s'            => \$species,
    "db=s"            => \$path_to_DB,
    "GO_file=s"       => \$GO_file,
    "plot"            => \$plot,
    "pval=s"          => \$pval_cutoff,
    "pval_type=s"     => \$meth,
    "max_plot=s"      => \$num_cutoff,
    "plot_only"       => \$plot_only,
    "open_plot"       => \$plot_open,
    "install_topGO"   => \$install_topGO,
    "filt_enrich=s"   => \$sort_enrich,
    "min_genes_req=s" => \$min_genes,
    "max_genes_req=s" => \$max_genes,
    "allow_inferred"  => \$infer,
);

#------------------------------------------------------------------------------
# Validate arguments
#------------------------------------------------------------------------------
my $EXIT_STATUS = 0;
if (defined $input && ((defined $species && defined $path_to_DB) || (defined $GO_file))) {
    print "Starting...\n\n";
}
else {
    $EXIT_STATUS = 1;
}

if ($helpFlag or $EXIT_STATUS) {
    die "
usage: ChopGO.pl -i <INPUT_list> (-sp <SPECIES_CODE> -db </path/to/DB> OR --GO_file Custom_GO_file) [options]

compulsory:
    -i               Input list of genes (compulsory). Sep by \'\\n\'. Optional second column for subdivisions of list.  

    -sp              Species must be (Hsa,Bla,Mmu code). Check availability in DB.
OR
    --GO_file        Custom file with GeneID\tGO_term associations

optional:
    -bg              Background must be a list of genes sep by \\n.
    -db              Database, give path to DB of GO annotations. Default (./VASTDB)
    -install_topGO   Install topGO in R.

plotting options:

    -plot            Plot (BOOLEAN, optional), will create a histogram/s. Default=0=(OFF).
    -pval            Choose pval cutoff for histogram plots. Default=0.05.
    -pval_type	     Choose correction (holm, hochberg, hommel, bonferroni, BH, BY, fdr or none). Default=none.
    -filt_enrich     Filter by fold enrichment. Default=0.
    -max_plot        Max number of results to plot for each histogram (Bp,MF,CC). Default=10.
    -min_genes_req   Choose minimum number of genes containing each GO term. Default=2.
    -max_genes_req   Choose maximum number of genes containing each GO term. Default=100000.
    -allow_inferred  Allow inferred parent (GO, not in original file) into final enrichment. Default=0=OFF.
    -plot_only       Plot Only (BOOLEAN, optional), don't do GO enrichments. Default=0=OFF.
    -open_plot       Open plots once made (BOOLEAN, optional). Default=0=OFF.

Pre-requisites:
    R 

*** Questions: 
    Chris Wyatt (cw13722@gmail.com)

";
}

#------------------------------------------------------------------------------
# Resolve the GeneID -> GO mapping file
#------------------------------------------------------------------------------
my $path_to_GOs;
if ($infer) {
    # inferred file path supplied by caller; nothing to resolve here
}
elsif ($GO_file) {
    $path_to_GOs = $GO_file;
    $infer       = $path_to_GOs;
}
else {
    $path_to_GOs = "$path_to_DB\/$species\/FILES/GO_FILE_$species";
    $infer       = $path_to_GOs;
}

if (-e $path_to_GOs) {
    # good for you
}
else {
    print "

GO file does not exist:\n
= $path_to_GOs\n
This may be because there is no GO file in the DB you are pointing at.\n
If this is the case, you might do the following: \n
1,  Go to https://github.com/vastgroup/vast-tools.
2,  Download the GO_file of the species of interest.

Alternatively, you may download your own annotation, e.g.:

1,  Go to http://www.ensembl.org/index.html.\n
2,  Go to BioMart.\n
3,  Choose database Ensembl Genes (latest version).\n
4,  Choose species.\n 
5,  Click Attributes.\n
6,  Under Gene, select \"Gene ID\" only.\n
7,  Under External, select : \"GO Term Accession\" and \"GO Term Name\".\n
8,  Click results, Then Go (to save the file to your computer).\n
9,  Finally, make sure file is tab separated, and it goes GENE ID (e.g. ENG000..), GO (e.g.GO:000001), GO name (e.g. Liver related). IN this order.\n
10, Save file as GO_FILE_species (e.g. GO_FILE_Hsa), in the DB folder under the same three letter species name.\n
	";
}

#==============================================================================
# MAIN: full enrichment (unless -plot_only was requested)
#==============================================================================
if ($plot_only == 0) {

    # Output tables generated by R, collected for the plotting step.
    my @ALL_made_files;

    # Multiple-testing correction methods written into the R script.
    my @methods = ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none");

    # Report the chosen GO source back to the user.
    print "GO = $species GO\nDb = $path_to_DB\n\n" if !defined $GO_file;
    print "GO file = $GO_file\n\n"                 if  defined $GO_file;

    #--------------------------------------------------------------------------
    # BACKGROUND
    #--------------------------------------------------------------------------
    my $outfile = "BACKGROUND\.forR";
    my %Background_hash;

    if ($background) {
        # User-supplied background: read the gene list, then keep only those
        # genes from the GO mapping.
        open(my $IN_b, "<", $background) or die "Could not open $background \n";
        while (my $line = <$IN_b>) {
            $line =~ s/\r//g;
            chomp $line;

            # Rename weird NCBI id rna- prefix
            if ($line =~ m/rna-/) {
                my @sp1 = split(/\-/, $line);
                $line = $sp1[1];
            }

            # Rename weird transcript ids with ':', usually transcript:ENSGMT0000012,
            # we want just ENSGMT0000012
            if ($line =~ m/\:/) {
                my @sp1 = split(/\:/, $line);
                $line = $sp1[1];
            }

            my @input_here = split("\t", $line);
            if (scalar @input_here > 1.5) {
                $Background_hash{$input_here[0]} = "HIT";
            }
            else {
                $Background_hash{$line} = "HIT";
            }
        }
        my $n_back = keys %Background_hash;

        # Build the gene2GO background for R (restricted to the -bg genes).
        my $Table_b1;
        if ($GO_file) {
            $Table_b1 = $GO_file;
        }
        else {
            $Table_b1 = "$path_to_DB\/$species\/FILES/GO_FILE_$species";
        }
        open(my $IN,        "<", $Table_b1) or die "Could not open $Table_b1 \n";
        open(my $outhandle, ">", $outfile)  or die "Could not open $outfile \n";
        my %Gene_Go_Hash;
        my %tot_genes;
        while (my $line = <$IN>) {
            $line =~ s/\r//g;
            chomp $line;
            my @linesplit = split("\t", $line);
            my $gene = $linesplit[0];
            $tot_genes{$gene} = "HIT";
            if ($Background_hash{$gene}) {
                my $GO = $linesplit[1];
                if ($GO) {
                    if ($Gene_Go_Hash{$gene}) {
                        my $old = $Gene_Go_Hash{$gene};
                        $Gene_Go_Hash{$gene} = "$old\",\"$GO";
                    }
                    else {
                        $Gene_Go_Hash{$gene} = $GO;
                    }
                }
            }
        }
        my $n_all = keys %tot_genes;
        print "BACKGROUND\nOf $n_all total Genes with GO annotation, $n_back were chosen for the background (based on user -b list)\n\n";
        print "variable^^^\n";
        print $outhandle "Chop.gene2GO<- list()\n";
        foreach my $key (keys %Gene_Go_Hash) {
            print $outhandle "Chop.gene2GO\$$key <- c(\"$Gene_Go_Hash{$key}\")\n";
        }
    }
    else {
        # No background chosen: use all annotated genes as the background.
        my $Table_b2;
        if ($GO_file) {
            $Table_b2 = $GO_file;
        }
        else {
            $Table_b2 = "$path_to_DB\/$species\/FILES/GO_FILE_$species";
        }
        open(my $IN,        "<", $Table_b2) or die "Could not open $Table_b2 \n";
        open(my $outhandle, ">", $outfile)  or die "Could not open $outfile \n";
        my %Gene_Go_Hash;
        while (my $line = <$IN>) {
            $line =~ s/\r//g;
            chomp $line;
            my @linesplit = split("\t", $line);
            my $gene = $linesplit[0];
            $gene =~ s/-/_/g;
            $gene =~ s/:/_/g;
            my $GO = $linesplit[1];
            if ($GO) {
                if ($Gene_Go_Hash{$gene}) {
                    my $old = $Gene_Go_Hash{$gene};
                    $Gene_Go_Hash{$gene} = "$old\",\"$GO";
                }
                else {
                    $Gene_Go_Hash{$gene} = $GO;
                }
            }
        }

        print $outhandle "Chop.gene2GO<- list()\n";
        foreach my $key (keys %Gene_Go_Hash) {
            print $outhandle "Chop.gene2GO\$$key <- c(\"$Gene_Go_Hash{$key}\")\n";
        }
    }

    #--------------------------------------------------------------------------
    # QUERY INPUT
    #--------------------------------------------------------------------------
    my $Table_i = $input;
    open(my $IN, "<", $Table_i) or die "Could not open $Table_i \n";

    # Strip any directory component from the file name used downstream in R.
    if ($Table_i =~ m/\//) {
        my @split = split("\/", $Table_i);
        $Table_i = $split[-1];
        $Table_i =~ s/\-/\_/g;
    }

    # R cannot use "-" in object names, so sanitise the name.
    if ($Table_i =~ m/\-/) {
        print "Name of input file contains \"-\"\'s , replacing with \"\_\"\'s, for processing in R\n";
        $Table_i =~ s/\-/\_/g;
    }

    my $outfile2 = "$Table_i\.converted";
    my $outfile3 = "$Table_i\.R_GO_subs";

    open(my $outhandle2, ">", $outfile2) or die "Could not open $outfile2 \n";
    open(my $outhandle3, ">", $outfile3) or die "Could not open $outfile3 \n";

    print "\nStep 2\n\n";

    # Group genes by the optional second column (or by file name if absent).
    my %Gene_Go_Hash;
    while (my $line = <$IN>) {
        $line =~ s/\r//g;
        $line =~ s/\"//g;
        chomp $line;
        my @linesplit = split("\t", $line);
        my $gene = $linesplit[0];
        $gene =~ s/-/_/g;
        $gene =~ s/:/_/g;

        if ($linesplit[1]) {
            my $GROUP = $linesplit[1];
            if ($Gene_Go_Hash{$GROUP}) {
                my $old = $Gene_Go_Hash{$GROUP};
                $Gene_Go_Hash{$GROUP} = "$old\",\"$gene";
            }
            else {
                $Gene_Go_Hash{$GROUP} = $gene;
            }
        }
        else {
            my $GROUP = "$Table_i";
            if ($Gene_Go_Hash{$GROUP}) {
                my $old = $Gene_Go_Hash{$GROUP};
                $Gene_Go_Hash{$GROUP} = "$old\",\"$gene";
            }
            else {
                $Gene_Go_Hash{$GROUP} = $gene;
            }
        }
    }

    #--------------------------------------------------------------------------
    # BUILD R COMPARISONS (one enrichment per group, across BP / MF / CC)
    #--------------------------------------------------------------------------
    print $outhandle2 "Chop.WGCNA2Gene<- list()\n";

    foreach my $key (keys %Gene_Go_Hash) {

        print $outhandle2 "Chop.WGCNA2Gene\$$key <- c(\"$Gene_Go_Hash{$key}\")\n";

        # Select genes for this group and flag them against the background.
        print $outhandle3 "selGenes<-Chop.WGCNA2Gene\$", $key, "\n";
        print $outhandle3 "inGenes <- factor(as.integer(names(Chop.gene2GO) %in% selGenes))\n";
        print $outhandle3 "names(inGenes) <- names(Chop.gene2GO)\n";

        # BP
        print $outhandle3 "GOdata <- new(\"topGOdata\", ontology=\"BP\", allGenes=inGenes, annot=annFUN.gene2GO, gene2GO=Chop.gene2GO)\n";
        print $outhandle3 "resultFisher <- runTest(GOdata, algorithm = \"classic\", statistic = \"fisher\")\n";
        print $outhandle3 "allRes_BP <- GenTable(GOdata, classicFisher = resultFisher,orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 50)\n";
        foreach my $meths (@methods) {
            print $outhandle3 "allRes_BP\$", $meths, "<-p.adjust(allRes_BP\$classicFisher, method = \"", $meths, "\")\n";
        }
        print $outhandle3 "allRes_BP\$FoldChange<-allRes_BP\$Significant/allRes_BP\$Expected\n";
        print $outhandle3 "allRes_BP\$ontology<-\"BP\"\n";

        # MF
        print $outhandle3 "GOdata <- new(\"topGOdata\", ontology=\"MF\", allGenes=inGenes, annot=annFUN.gene2GO, gene2GO=Chop.gene2GO)\n";
        print $outhandle3 "resultFisher <- runTest(GOdata, algorithm = \"classic\", statistic = \"fisher\")\n";
        print $outhandle3 "allRes_MF <- GenTable(GOdata, classicFisher = resultFisher,orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 50)\n";
        foreach my $meths (@methods) {
            print $outhandle3 "allRes_MF\$", $meths, "<-p.adjust(allRes_MF\$classicFisher, method = \"", $meths, "\")\n";
        }
        print $outhandle3 "allRes_MF\$FoldChange<-allRes_MF\$Significant/allRes_MF\$Expected\n";
        print $outhandle3 "allRes_MF\$ontology<-\"MF\"\n";

        # CC
        print $outhandle3 "GOdata <- new(\"topGOdata\", ontology=\"CC\", allGenes=inGenes, annot=annFUN.gene2GO, gene2GO=Chop.gene2GO)\n";
        print $outhandle3 "resultFisher <- runTest(GOdata, algorithm = \"classic\", statistic = \"fisher\")\n";
        print $outhandle3 "allRes_CC <- GenTable(GOdata, classicFisher = resultFisher,orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 50)\n";
        foreach my $meths (@methods) {
            print $outhandle3 "allRes_CC\$", $meths, "<-p.adjust(allRes_CC\$classicFisher, method = \"", $meths, "\")\n";
        }
        print $outhandle3 "allRes_CC\$FoldChange<-allRes_CC\$Significant/allRes_CC\$Expected\n";
        print $outhandle3 "allRes_CC\$ontology<-\"CC\"\n";

        # Merge ontologies, filter by p-value / fold enrichment, and write out.
        print $outhandle3 "ALL_res_", $key, "<-rbind(allRes_BP, allRes_MF, allRes_CC)\n";
        print $outhandle3 "x.sub <- subset(ALL_res_", $key, ", none < ", $pval_cutoff, ")\n";
        print $outhandle3 "f.sub<-x.sub[sort.list(x.sub\$none),]\n";
        print $outhandle3 "e.sub <- subset(f.sub, FoldChange > ", $sort_enrich, ")\n";
        print $outhandle3 "write.table(e.sub, \"", $key, "_TopGo_results_ALL.tab\", sep=\"\\t\", quote=FALSE, eol=\"\\n\", row.names=F)   \n";
        my $out_name = "$key\_TopGo_results_ALL.tab";
        push(@ALL_made_files, $out_name);
    }

    print "Finished compiling scripts, now running R... (can take ~40 seconds for each query)\n\n";

    #--------------------------------------------------------------------------
    # RUN R
    #--------------------------------------------------------------------------
    my $outfile4 = "R_CODE_TO_RUN";
    open(my $outhandle4, ">", $outfile4) or die "Could not open $outfile4 \n";
    if (defined $install_topGO) {
        print $outhandle4 "source(\"http://bioconductor.org/biocLite.R\")\n";
        print $outhandle4 "biocLite(\"topGO\")\n";
    }

    print $outhandle4 "suppressWarnings(suppressMessages(library (topGO)))\n";
    print $outhandle4 "source (\"$outfile\")\nsource (\"$outfile2\")\nsource (\"$outfile3\")\n";
    print $outhandle4 "save.image(\"Image.Rdata\")\nquit(save = \"no\")\n";

    `R --vanilla < R_CODE_TO_RUN >output.ofthis.test`;

    print "Start Plotting\n\n";

    #--------------------------------------------------------------------------
    # PLOT HISTOGRAMS
    #--------------------------------------------------------------------------
    if (defined $plot) {
        my @files = @ALL_made_files;
        foreach my $results (@files) {
            chomp $results;
            my $in_here = "$binPath\/MakeHist_fromChartGO-10.pl";
            if ($in_here =~ m/ /) {
                $in_here =~ s/ /\\ /g;
            }
            if ($in_here =~ m/\(/) {
                $in_here =~ s/\(/\\(/g;
            }
            if ($in_here =~ m/\)/) {
                $in_here =~ s/\)/\\)/g;
            }
            print "HERE: $in_here $results $pval_cutoff $num_cutoff $plot_open $meth $sort_enrich $min_genes $max_genes $infer\n";
            `$in_here $results $pval_cutoff $num_cutoff $plot_open $meth $sort_enrich $min_genes $max_genes $infer`;
        }
    }

}
#==============================================================================
# PLOT ONLY: skip enrichment and just (re)generate the histograms
#==============================================================================
else {
    my $Table_i = $input;
    open(my $IN, "<", $Table_i) or die "Could not open $Table_i \n";
    my %Gene_Go_Hash;
    while (my $line = <$IN>) {
        $line =~ s/\r//g;
        $line =~ s/\"//g;
        chomp $line;
        my @linesplit = split("\t", $line);
        my $gene = $linesplit[0];
        if ($linesplit[1]) {
            my $GROUP = $linesplit[1];
            if ($Gene_Go_Hash{$GROUP}) {
                my $old = $Gene_Go_Hash{$GROUP};
                $Gene_Go_Hash{$GROUP} = "$old\",\"$gene";
            }
            else {
                $Gene_Go_Hash{$GROUP} = $gene;
            }
        }
        else {
            my $GROUP = "$Table_i";
            if ($Gene_Go_Hash{$GROUP}) {
                my $old = $Gene_Go_Hash{$GROUP};
                $Gene_Go_Hash{$GROUP} = "$old\",\"$gene";
            }
            else {
                $Gene_Go_Hash{$GROUP} = $gene;
            }
        }
    }

    my @ALL_made_files;
    foreach my $key (keys %Gene_Go_Hash) {
        my $out_name = "$key\_TopGo_results_ALL.tab";
        push(@ALL_made_files, $out_name);
    }

    foreach my $results (@ALL_made_files) {
        chomp $results;
        my $in_here = "$binPath\/MakeHist_fromChartGO-10.pl";
        if ($in_here =~ m/ /) {
            $in_here =~ s/ /\\ /g;
        }
        if ($in_here =~ m/\(/) {
            $in_here =~ s/\(/\\(/g;
        }
        if ($in_here =~ m/\)/) {
            $in_here =~ s/\)/\\)/g;
        }
        print "$in_here $results $pval_cutoff $num_cutoff $plot_open $meth $sort_enrich $min_genes $max_genes $infer\n";
        `perl $in_here $results $pval_cutoff $num_cutoff $plot_open $meth $sort_enrich $min_genes $max_genes $infer`;
    }
}

print "Script completed\n\n";

# Tidy up temporary file naming.
$input =~ s/\-/\_/g;