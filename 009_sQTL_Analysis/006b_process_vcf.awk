BEGIN{OFS="\t";FS="\t"}
{
    l=$1"\t"$2"\t"$3;
    for(i=10;i<=NF;i++){
    if($1=="#CHROM"){
        g=$i
    } else {
        g="NA"
        if($i~"0/0:"){g=0}
        if($i~"1/0:"){g=1}
        if($i~"0/1:"){g=1}
        if($i~"1/1:"){g=2}
    }
    l=l"\t"g
    }
    print l;
}