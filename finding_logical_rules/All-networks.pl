#!/usr/bin/perl;

###
### Developed by: Juan A. Arias Del Angel (Universidad Nacional Autónoma de México)
### Advisors: Mariana Benitez, Ana E. Escalante and Leon Patricio Martinez Castilla.
### Reference: Arias Del Angel, Martinez-Castilla, Escalante and Benitez (in preparation).
### 
### Description: Taking a first list of "confirmed" directed interactions and a second one of "putative" non-directed interactons,
### the script performs the following steps:
### 1) Generate all the possible subsets of "putative" interactions (i.e., containing from 0 to all the interactions).
### 2) For each subset of putative interactions, generate directed interactions.
### 3) Merge the "putative" directed interactions with the "confirmed" interactions.
### 4) Filter the resulting networks by the criteria that each node in the network has non zero in- and out-degree.
### 5) For each resulting network from the previous step, create dynamic boolean networks by assigning boolean functions.
### 6) Print the results.
### NOTE: Boolean function assignment is according the principles described in
### Azpeitia et al. (2013) Finding missing interactions of the Arabidopsis thaliana root stem cell niche gene regulatory network. 
### Possible boolean functions are listed in the files "Rules-X-Regulators.txt".
### Input: this script do not take any argument as input but some files are required to be in the same directory.
### Required files: 
### a) Rules-1-Regulators.txt
### b) Rules-2-Regulators.txt
### c) Rules-3-Regulators.txt
### d) Rules-4-Regulators.txt
### e) showAttractors.r
### Output file: the script print to screen the results. It is recommended to sent the output to an output file using the proper commands.
### For each of the tested network, the output consists of two rows with the following format:
### ROW 1 => RED N NODE 1 = REG1 REG2 ... REGN ( truth table 1 ) ; NODE 2 = REG1 REG2 ... REGN ( truth table 1 ) ; ...
### ROW 2 => This networks has K attractors: Attr1//Attr2//Attr3.1-Attr3.2-...-Attr3.n//...
### 
### In ROW 2 attractors are codified in decimal notation. Attractors are separated by '//' and states of a single cyclic attractor 
### are separated by '-'. 
### Note that the code is building by using subroutines and so, the reading of the code is not be lineal. 

use strict;


###########################################################################################
###                                                                                     ###
###                                     MAIN SCRIPT                                     ###
###                                                                                     ###
###########################################################################################


### List of "confirmed" interactions. 

my @current_interactions = ("STARV --> STARV", "STARV --> PEP", "PEP --> PKTD1", "PKTC2 --> PSKA5", "PSKA5 --> MRPC2", "MRPC2 --> PKTA2", "DEVRS --> PKTA2", "STARV --> PKTD9", "MRPC2 --> FRUA", "MRPC2 --> DEVRS", "FRUA --> DEVRS", "DEVRS --> DEVRS", "FRUA --> FRUA", "DEVRS --> FRUA");

### List of "putative" interactins. 
my @putative_interactions = ('MKAPA --> PKTA2', 'MKAPA --> PKTA4', 'MKAPA --> PKTC2', 'MKAPA --> PKTD1', 'MKAPA --> PKTD9', 'MKAPB --> PKTA4', 'MKAPB --> PKTC2', 'MKAPB --> PKTD1', 'MKAPC --> PKTC2', 'MKAPC --> PKTD1');


my $net_id = 0;

&create_subset(@putative_interactions);

########################################################################################################################
########################################################################################################################


### Generate all the possible subsets of "putative" interactions.
sub create_subset {
    my @interactions = @_;
    my $N = scalar(@interactions);
    
    for(my $k = 0; $k < 2**$N; $k++){
        my @bits = ((0) x $N);
        my $bin = sprintf ("%b",$k);
        my @bin = split(//, $bin);
        my @pos = ((0) x ($N - scalar(@bin)), @bin);
        
        my @subset = ();
        for(my $i = 0; $i < $N; $i++){
            if($pos[$i] == 1){
                push(@subset, $interactions[$i]);
            }
        }
        &assign_direction(@subset);
        
    }
}

### For each subset, assign direction to every interaction.
sub assign_direction {
    my @interactions = @_;
    my $N = scalar(@interactions);
    my @elements = ();
    my $interaction = ();
    
    for(my $k = 0; $k < 2**$N; $k++){
        my @bits = ((0) x $N);
        my $bin = sprintf ("%b",$k);
        my @bin = split(//, $bin);
        my @pos = ((0) x ($N - scalar(@bin)), @bin);
        
        my @set = ();
        for(my $i = 0; $i < $N; $i++){
            if($pos[$i] == 1){
                @elements = split(' --> ', $interactions[$i]);
                $interaction = join(' --> ', $elements[1], $elements[0]);
                push(@set, $interaction);
            }
            else {
                push(@set, $interactions[$i])
            }
        }
        push(@set, @current_interactions);
        &filter_by_io(@set);
        
    }
}

### Filter directed networks by in- and -out degree. 
sub filter_by_io {
    my @net = @_;
    my %input = ();
    my %output = ();
    my $interaction = ();
    my @elements = ();
    my $in_nodes = ();
    my $out_nodes = ();
    
    for($interaction = 0; $interaction < scalar(@net); $interaction++){
        @elements = split(' --> ', $net[$interaction]);
        $input{$elements[0]} += 1;
        $output{$elements[1]} += 1;
        
    }
    
    $in_nodes = join(' ', sort keys %input);
    $out_nodes = join(' ', sort keys %output);
    
    if($in_nodes eq $out_nodes){
       $net_id++;
       if($net_id > 0){
         my $r_file = constructBoolNet(@net);
      	 open(my $R, '>', 'runR_boolnet.r');
         print $R $r_file;
       		
	 # Run a single Boolean model. 
       	 system("Rscript runR_boolnet.r $net_id");
	}
    }
    
    
    
}

### Generate Boolean models for each of the directed networks.
sub constructBoolNet{
    my @net = @_;
    
    my %in = ();
    my $rules = ();
    
    foreach my $interaction (@net){
        my @nodes = split(" \-\-\> ", $interaction);
        if($nodes[1] eq 'PKTA4' || $nodes[1] eq 'PKTC2' || $nodes[1] eq 'PKTD1' || $nodes[1] eq 'MKAPA' || $nodes[1] eq 'MKAPB' || $nodes[1] eq 'MKAPC'){
        push(@{$in{$nodes[1]}}, $nodes[0]);
        }
    }
    
    foreach my $node (sort keys %in){
        #print "$node, @{$in{$node}}\n";
    }
     my %pep_rules = (1 => "pep.rules = subset(one.r, one.r[,3] == \'-\')\n\n", 2 => "pep.rules = subset(two.r, two.r[,6] == \'-\')\n\n", 3 => "pep.rules = subset(three.r, three.r[,11] == \'-\')\n\n");
    
    
    my $file = "library(BoolNet)\n";
    $file .= "source\(\"showAttractors.r\"\)\n\n";
    $file .= "args = commandArgs(TRUE)\n";
    $file .= "id = args\[1\]\n\n";
    $file .= "one.r   = read.csv\(\"Rules-1-Regulators.txt\", header = F, sep =\" \"\)\n";
    $file .= "two.r   = read.csv\(\"Rules-2-Regulators.txt\", header = F, sep =\" \"\)\n";
    $file .= "three.r = read.csv\(\"Rules-3-Regulators.txt\", header = F, sep =\" \"\)\n";
    $file .= "four.r  = read.csv\(\"Rules-4-Regulators.txt\", header = F, sep =\" \"\)\n";
    
    $file .= $pep_rules{scalar(@{$in{'PKTD1'}})};
    
    
    $file .= "sink\(\"testNet.bn\"\)\n";
    $file .= "cat\(\"targets, factors\\n\"\)\n";
    
    foreach my $node (sort keys %in){
        $file .= "cat\(\"" . $node . ", " . join(" \& ", @{$in{$node}}) . "\\n\"\)\n";
    }
    
    $file .= "cat\(\"PKTA2, MRPC2\\n\"\)\n";
    $file .= "cat\(\"PKTD9, !STARV\\n\"\)\n";
    $file .= "cat\(\"PSKA5, PKTC2\\n\"\)\n";
    $file .= "cat\(\"PEP, !STARV\\n\"\)\n";
    $file .= "cat\(\"STARV, STARV\\n\"\)\n";
    $file .= "cat\(\"MRPC2, !PSKA5 \| MRPC2\\n\"\)\n";
    $file .= "sink\(\)\n\n";
    
    $file .= "net = loadNetwork\(\"testNet.bn\"\)\n";
    
    my $nodes = join(' ', sort keys %in);
    
    if($nodes =~ /MKAPA/){
        if(scalar(@{$in{'MKAPA'}}) == 1){
        $file .= "for\(mkapA in 1\:nrow\(one\.r\)\)\{\n";
        $file .= "\tnet\$interactions\$MKAPA\$func = one\.r\[mkapA, 1\:2\]\n";
        }
        if(scalar(@{$in{'MKAPA'}}) == 2){
            $file .= "for\(mkapA in 1\:nrow\(two\.r\)\)\{\n";
            $file .= "\tnet\$interactions\$MKAPA\$func = two\.r\[mkapA, 1\:4\]\n";
        }
        if(scalar(@{$in{'MKAPA'}}) == 3){
            $file .= "for\(mkapA in 1\:nrow\(three\.r\)\)\{\n";
            $file .= "\tnet\$interactions\$MKAPA\$func = three\.r\[mkapA, 1\:8\]\n";
        }
        if(scalar(@{$in{'MKAPA'}}) == 4){
            $file .= "for\(mkapA in 1\:nrow\(four\.r\)\)\{\n";
            $file .= "\tnet\$interactions\$MKAPA\$func = four\.r\[mkapA, 1\:16\]\n";
        }
    }
    
    if($nodes =~ /MKAPB/){
        if(scalar(@{$in{'MKAPB'}}) == 1){
        $file .= "\tfor\(mkapB in 1\:nrow\(one\.r\)\)\{\n";
        $file .= "\t\tnet\$interactions\$MKAPB\$func = one\.r\[mkapB, 1\:2\]\n";
        }
        if(scalar(@{$in{'MKAPB'}}) == 2){
            $file .= "\tfor\(mkapB in 1\:nrow\(two\.r\)\)\{\n";
            $file .= "\t\tnet\$interactions\$MKAPB\$func = two\.r\[mkapB, 1\:4\]\n";
        }
        if(scalar(@{$in{'MKAPB'}}) == 3){
            $file .= "\tfor\(mkapB in 1\:nrow\(three\.r\)\)\{\n";
            $file .= "\t\tnet\$interactions\$MKAPB\$func = three\.r\[mkapB, 1\:8\]\n";
        }
    }
    
    if($nodes =~ /PKTC2/){
        if(scalar(@{$in{'PKTC2'}}) == 1){
            $file .= "\t\tfor\(pktc2 in 1\:nrow\(one\.r\)\)\{\n";
            $file .= "\t\t\tnet\$interactions\$PKTC2\$func = one\.r\[pktc2, 1\:2\]\n";
        }
        if(scalar(@{$in{'PKTC2'}}) == 2){
            $file .= "\t\tfor\(pktc2 in 1\:nrow\(two\.r\)\)\{\n";
            $file .= "\t\t\tnet\$interactions\$PKTC2\$func = two\.r\[pktc2, 1\:4\]\n";
        }
        if(scalar(@{$in{'PKTC2'}}) == 3){
            $file .= "\t\tfor\(pktc2 in 1\:nrow\(three\.r\)\)\{\n";
            $file .= "\t\t\tnet\$interactions\$PKTC2\$func = three\.r\[pktc2, 1\:8\]\n";
        }
    }
    
    if($nodes =~ /PKTD1/){
        if(scalar(@{$in{'PKTD1'}}) == 1){
            $file .= "\t\t\tfor\(pktd1 in 1\:nrow\(pep\.rules\)\)\{\n";
            $file .= "\t\t\t\tnet\$interactions\$PKTD1\$func = pep\.rules\[pktd1, 1\:2\]\n";
        }
        if(scalar(@{$in{'PKTD1'}}) == 2){
            $file .= "\t\t\tfor\(pktd1 in 1\:nrow\(pep\.rules\)\)\{\n";
            $file .= "\t\t\t\tnet\$interactions\$PKTD1\$func = pep\.rules\[pktd1, 1\:4\]\n";
        }
        if(scalar(@{$in{'PKTD1'}}) == 3){
            $file .= "\t\t\tfor\(pktd1 in 1\:nrow\(pep\.rules\)\)\{\n";
            $file .= "\t\t\t\tnet\$interactions\$PKTD1\$func = pep\.rules\[pktd1, 1\:8\]\n";
        }
    }
    
    $file .= "\t\t\t\tattr = getAttractors\(net\)\n";
    $file .= "\t\t\t\tprintNet\(net, id\)\n";
    $file .= "\t\t\t\tshowAttractors\(attr\)\n";
    if($nodes =~ /PKTD1/){
        $file .= "\t\t\t\}\n";
    }
    if($nodes =~ /PKTC2/){
        $file .= "\t\t\}\n";
    }
    if($nodes =~ /MKAPB/){
        $file .= "\t\}\n";
    }
    if($nodes =~ /MKAPA/){
        $file .= "\}\n";
    }
    return($file);

}
