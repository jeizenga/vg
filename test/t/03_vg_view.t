#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 20

is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view -d - | wc -l) 505 "view produces the expected number of lines of dot output"
is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view -g - | wc -l) 503 "view produces the expected number of lines of GFA output"
# This one may throw warnings related to dangling edges because we just take an arbitrary subset of the GFA
is $(vg construct -r small/x.fa -v small/x.vcf.gz | vg view - | head | vg view -Fv - | vg view - | wc -l) 10 "view converts back and forth between GFA and vg format"

is $(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -a - | wc -l) $(samtools view -u minigiab/NA12878.chr22.tiny.bam | samtools view - | wc -l) "view can convert BAM to GAM"

is "$(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -aj - | jq -c --sort-keys . | sort | md5sum)" "$(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -aj - | vg view -JGa - | vg view -aj - | jq -c --sort-keys . | sort | md5sum)" "view can round-trip JSON and GAM"

# We need to run through GFA because vg construct doesn't necessarily chunk the
# graph the way vg view wants to.
vg construct -r small/x.fa -v small/x.vcf.gz | vg view -g - | vg view -Fv - >x.vg
vg view -j x.vg | jq . | vg view -Jv - | diff x.vg -
is $? 0 "view can reconstruct a VG graph from JSON"

vg view -v x.vg | cmp -s - x.vg
is $? 0 "view can pass through VG"

rm -f x.vg

is $(samtools view -u minigiab/NA12878.chr22.tiny.bam | vg view -bG - | vg view -a - | jq .sample_name | grep -v '^"1"$' | wc -l ) 0 "view parses sample names"

is $(vg view -f ./small/x.fa_1.fastq  ./small/x.fa_2.fastq | vg view -a - | wc -l) 2000 "view can handle fastq input"

is $(vg view -Jv ./cyclic/two_node.json | vg view -j - | jq ".edge | length") 4 "view can translate graphs with 2-node cycles"

is $(vg view -g ./cyclic/all.vg | tr '\t' ' ' | grep "4 + 4 -" | wc -l) 1 "view outputs properly oriented GFA"

is $(vg view -d ./cyclic/all.vg | wc -l) 23 "view produces the expected number of lines of dot output from a cyclic graph"

# We need to make a single-chunk graph
vg construct -r small/x.fa -v small/x.vcf.gz | vg view -v - >x.vg
is $(cat x.vg x.vg x.vg x.vg | vg view -c - | wc -l) 4 "streaming JSON output produces the expected number of chunks"

is "$(cat x.vg x.vg | vg view -vD - 2>&1 > /dev/null | wc -l)" 0 "duplicate warnings can be suppressed"

rm x.vg

is "$(vg view -Fv overlaps/two_snvs_assembly1.gfa | vg stats -l - | cut -f2)" "315" "gfa graphs are imported pre-bluntified"

is "$(vg view -Fv overlaps/two_snvs_assembly1.gfa | vg mod --bluntify - | vg stats -l - | cut -f2)" "315" "bluntifying has no effect"

is "$(vg view -Fv overlaps/two_snvs_assembly4.gfa | vg stats -l - | cut -f2)" "335" "a more complex GFA can be imported"

vg view -Fv overlaps/incorrect_overlap.gfa >/dev/null 2>errors.txt
is "$?" "1" "GFA import rejects a GFA file with an overlap that goes beyond its sequences"
is "$(cat errors.txt | wc -l)" "1" "GFA import produces a concise error message in that case"

rm -f errors.txt

vg view -Fv overlaps/corrected_overlap.gfa >/dev/null
is "$?" "0" "GFA import accepts that file when the offending overlap length is fixed"

